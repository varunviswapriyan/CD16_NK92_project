"""
replicate_ci_gpr.py - Replicate Caldara & Iacoviello (2022, AER) GPR / GPT / GPA
through the same ProQuest TDM Studio search API used by discover_facets.py
(POST /api/cm/document/search, reading `docsFound`).

Two corpora, as in C&I:

  hist    Historical index (C&I: GPRH, GPRHT, GPRHA)
          New York Times, Washington Post, Chicago Tribune
          (ProQuest Historical Newspapers databases), 1900 onward,
          normalized to mean 100 over 1900-2019.

  recent  Recent index (C&I: GPR, GPRT, GPRA)
          C&I's full 10-paper panel: NYT, WaPo, Chicago Tribune, WSJ, Guardian,
          USA Today, LA Times, Daily Telegraph, Globe and Mail (by ProQuest
          publication ID) + the Financial Times. FT has no known publication ID
          in this subscription, so it is matched by title (PUB("Financial Times"))
          and/or any IDs you put in FT_PUBIDS; `ft recent` shows what ProQuest
          has for it. If FT returns nothing the index silently reduces to the
          9 papers; `probe` and `run` say so loudly.
          1985 onward, normalized to mean 100 over 1985-2019.

Method (matches C&I's published series exactly; checked against their
benchmark file: GPR = SHARE / mean(SHARE over ref period) * 100, and GPT/GPA
are each normalized to their own mean 100):

    share_t  = (# articles in month t matching the category query)
             / (# articles in month t matching the stopword denominator query)
    index_t  = share_t / mean(share over reference period) * 100

Numerator and denominator are always counted on the SAME corpus with the SAME
filters (sourcetype = Newspapers + the DTYPE clause), monthly. Because the API
returns exact counts, there is no 2M-document export limit and no annual
denominator interpolation (open item 3 in the methodology doc is resolved).

Query: C&I Table 1 as adopted in the consolidated methodology ("Build A",
which gave hist r = 0.982 / 0.961 / 0.989): "nuclear weapon*" included,
terror* (not terroris*), "offensive" in cat. 7 only, no "wage*", full
exclusion list through "tax".

Usage
-----
    export TDM_TOKEN="<bearer token>"

    python3 replicate_ci_gpr.py probe  hist          # 1-minute sanity checks, run FIRST
    python3 replicate_ci_gpr.py probe  recent
    python3 replicate_ci_gpr.py ft     recent        # find Financial Times titles / coverage
    python3 replicate_ci_gpr.py search recent --pattern financial  # any pub whose name has "financial"
    python3 replicate_ci_gpr.py titles recent        # list exact pubTitle names (debug aid)
    python3 replicate_ci_gpr.py run    hist recent   # build both indices (resumable)
    python3 replicate_ci_gpr.py validate hist recent --bench data_gpr_export.xls

`run` validates automatically when --bench is given. The benchmark can be
C&I's data_gpr_export.xls (needs xlrd) or the same sheet saved as .csv.

Outputs (under ./ci_replication/<corpus>/):
    monthly_counts.csv                       raw counts per month
    GPR_index.csv, GPT_index.csv, GPA_index.csv  (same schema as discover_facets.py)
    indices.png
    validation/comparison.csv, validation/report.txt, validation/overlay.png
"""

import argparse
import csv
import json
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date, datetime, timedelta

import requests

# ---------------------------------------------------------------------------
# API settings (same endpoint / payload shape as discover_facets.py)
# ---------------------------------------------------------------------------

TOKEN = os.environ.get("TDM_TOKEN", "<TOKEN>")
WORKBENCH_ID = os.environ.get("TDM_WORKBENCH_ID", "46334")
BASE_URL = "https://tdmstudio.proquest.com"

headers = {
    "Authorization": f"Bearer {TOKEN}",
    "Accept": "application/json, text/plain, */*",
    "Content-Type": "application/json",
}

REQUEST_PAUSE = 0.4
MAX_WORKERS = 8          # lower if you see repeated 429s
OUTPUT_ROOT = "ci_replication"

# ---------------------------------------------------------------------------
# Products (ProQuest databases). discover_facets.py used "globalnews".
# If `probe` reports 0 documents for a product, open the TDM workbench in a
# browser, run any search with that database ticked, and copy the exact
# {"moniker", "name"} pair from the /api/cm/document/search request payload
# (DevTools > Network). Paste it below.
# ---------------------------------------------------------------------------

HNP_NYT = {"moniker": "hnpnewyorktimes",
           "name": "ProQuest Historical Newspapers: The New York Times with Index"}
HNP_WAPO = {"moniker": "hnpwashingtonpost",
            "name": "ProQuest Historical Newspapers: The Washington Post"}
HNP_CHITRIB = {"moniker": "hnpchicagotribune",
               "name": "ProQuest Historical Newspapers: Chicago Tribune"}
GLOBALNEWS = {"moniker": "globalnews", "name": "Global Newsstream Collection"}

HIST_PRODUCTS = [HNP_NYT, HNP_WAPO, HNP_CHITRIB]
# Current-edition papers. Add the separate current NYT / WaPo / Chicago Tribune
# databases here too if `probe recent` shows those publication IDs returning 0
# under Global Newsstream alone.
RECENT_PRODUCTS = [GLOBALNEWS]

# ---------------------------------------------------------------------------
# Publication IDs (from Varun's methodology notes; Build B)
# ---------------------------------------------------------------------------

HIST_ARCHIVE_PUBIDS = {
    "New York Times":  ["55428"],
    "Washington Post": ["60400", "47009", "47011", "47012"],
    "Chicago Tribune": ["55409", "46851"],
}

RECENT_PUBIDS = {}  # unused — PUBID() does not resolve in Global Newsstream for these papers;
                    # all matching is now by exact publication title (confirmed via `search`).

# ---------------------------------------------------------------------------
# Recent-corpus publication titles — the EXACT names ProQuest returned for each
# paper via `search recent --pattern <word>`. Multiple titles per paper when
# ProQuest splits a paper into pre-/post-1997, Online, etc.: all editions of
# the same paper are OR'd together.
#
# Confirmed in this subscription (globalnews):
#   NYT:         "New York Times" (1985 onward, 596k articles total)
#   WaPo:        "The Washington Post" + "The Washington Post (pre-1997 Fulltext)"
#   ChiTrib:     "Chicago Tribune" + "Chicago Tribune (pre-1997 Fulltext)"
#   WSJ:         "Wall Street Journal" (US edition; Asian/Americas editions excluded by default
#                — C&I use the US WSJ only. Add them below if you want the broader panel.)
#   Guardian:    "The Guardian" + "The Guardian (pre-1997 Fulltext)" (UK Guardian; "Guardian"
#                alone also appears but is ambiguous — the many regional "...Guardian" papers
#                are correctly EXCLUDED by using PUB.EXACT instead of PUB.)
#   USA Today:   "USA TODAY" + "USA TODAY (pre-1997 Fulltext)" + "USA Today (Online)"
#   LA Times:    "Los Angeles Times" + "Los Angeles Times (pre-1997 Fulltext)"
#   Daily Tel.:  "The Daily Telegraph" + "The Sunday Telegraph" + "Telegraph.co.uk" +
#                "The Daily Telegraph (Online)"
#   Globe&Mail:  "The Globe and Mail"
#
# To add/remove a title, drop its line here and re-run `probe recent`.
# ---------------------------------------------------------------------------

RECENT_PAPER_TITLES = {
    "New York Times": [
        "New York Times",
    ],
    "Washington Post": [
        "The Washington Post",
        "The Washington Post (pre-1997 Fulltext)",
    ],
    "Chicago Tribune": [
        "Chicago Tribune",
        "Chicago Tribune (pre-1997 Fulltext)",
    ],
    "Wall Street Journal": [
        "Wall Street Journal",
    ],
    "The Guardian": [
        "The Guardian",
        "The Guardian (pre-1997 Fulltext)",
    ],
    "USA Today": [
        "USA TODAY",
        "USA TODAY (pre-1997 Fulltext)",
        "USA Today (Online)",
    ],
    "Los Angeles Times": [
        "Los Angeles Times",
        "Los Angeles Times (pre-1997 Fulltext)",
    ],
    "Daily Telegraph": [
        "The Daily Telegraph",
        "The Sunday Telegraph",
        "Telegraph.co.uk",
        "The Daily Telegraph (Online)",
    ],
    "Globe and Mail": [
        "The Globe and Mail",
    ],
}

# ---------------------------------------------------------------------------
# Financial Times — C&I's 10th paper. Confirmed via `search recent --pattern
# financial` and `--pattern FT` that FT is NOT in this subscription under any
# name (every match was Financial Express, Australian Financial Review, etc.,
# none of which are FT). Set INCLUDE_FT = True only after adding its database
# to FT_PRODUCTS and confirming via `ft recent`.
# ---------------------------------------------------------------------------
FT_PUBIDS = []
FT_TITLE_CLAUSES = ['PUB.EXACT("Financial Times")']  # tighter than PUB(...)
FT_PRODUCTS = []
INCLUDE_FT = False
REQUIRE_FT = False

RECENT_PAPERS = {
    name: {"pubids": [], "clauses": [f'PUB.EXACT("{t}")' for t in titles]}
    for name, titles in RECENT_PAPER_TITLES.items()
}
if INCLUDE_FT:
    RECENT_PAPERS["Financial Times"] = {"pubids": FT_PUBIDS, "clauses": FT_TITLE_CLAUSES}

# How a publication restriction is applied:
#   "query"    -> appends AND (PUBID(x) OR PUBID(y) ...) to every query (default)
#   "pubtitle" -> uses the pubTitle facet filter, as discover_facets.py does;
#                 fill RECENT_PUBTITLES with exact names from the `titles` command
#                 (include the FT title shown by `ft recent`)
PUB_FILTER_MODE = "query"
RECENT_PUBTITLES = []

# Historical corpus mode:
#   False -> ProQuest Historical Newspapers archives only (Build A, best match;
#            archives thin out after ~2008: WaPo stops ~2008, ChiTrib ~2015)
#   True  -> archives + current editions of the same 3 papers, by publication ID,
#            so the 3-paper series runs to the present (Build B approach)
HIST_EXTEND_WITH_CURRENT = False

# Also count the 8 categories separately (8 extra calls per month). Lets you
# compare against C&I's SHAREH_CAT_1..8 to see which category drives any gap.
CATEGORY_COUNTS = False

# ---------------------------------------------------------------------------
# Search dictionary - C&I Table 1, as adopted ("Build A")
# ---------------------------------------------------------------------------

DTYPE_PREFIX = (
    'DTYPE(article OR commentary OR editorial OR feature OR '
    '"front page article" OR "front page/cover story" OR news OR report OR review)'
)

_THREAT = ('risk* OR warn* OR fear* OR danger* OR threat* OR doubt* OR crisis OR '
           'troubl* OR disput* OR concern* OR tension* OR imminen* OR inevitable OR '
           'footing OR menace* OR brink OR scare OR peril*')
_WAR = ('war OR conflict OR hostilities OR revolution* OR insurrection OR uprising OR '
        'revolt OR coup OR geopolitical')
_TERROR = 'terror* OR guerrilla* OR hostage*'
_NUCLEAR = (
    '("nuclear war" OR "nuclear warfare" OR "nuclear warhead" OR "nuclear warheads" OR "nuclear wars") OR '
    '("atomic war" OR "atomic warfare" OR "atomic warheads" OR "atomic wars") OR '
    '("nuclear missile" OR "nuclear missiles") OR '
    '("nuclear bomb" OR "nuclear bombardment" OR "nuclear bomber" OR "nuclear bombers" OR '
    '"nuclear bombing" OR "nuclear bombs") OR '
    '("atomic bomb" OR "atomic bombing" OR "atomic bombings" OR "atomic bombs") OR '
    '"h-bomb*" OR ("hydrogen bomb" OR "hydrogen bombs") OR "nuclear test" OR "nuclear weapon*"'
)

CATEGORIES = {
    1: f'({_WAR}) NEAR/2 ({_THREAT})',
    2: '(peace OR truce OR armistice OR treaty OR parley) NEAR/2 '
       '(menace* OR reject* OR threat* OR peril* OR boycott* OR disrupt*)',
    3: '(military OR troops OR missile* OR "arms" OR weapon* OR bomb* OR warhead*) AND '
       '(buildup* OR "build-up*" OR blockad* OR sanction* OR embargo OR quarantine OR '
       'ultimatum OR mobiliz*)',
    4: f'({_NUCLEAR}) AND ({_THREAT})',
    5: f'({_TERROR}) NEAR/2 ({_THREAT})',
    6: f'({_WAR}) NEAR/2 (begin* OR begun OR began OR outbreak OR "broke out" OR breakout OR '
       'start* OR declar* OR proclamation OR launch*)',
    7: '(allie* OR enem* OR foe* OR army OR navy OR aerial OR troops OR rebels OR insurgen*) '
       'NEAR/2 (drive* OR shell* OR advance* OR offensive OR invasion OR invad* OR clash* OR '
       'attack* OR raid* OR launch* OR strike*)',
    8: f'({_TERROR}) NEAR/2 (act OR attack OR bomb* OR kill* OR strike* OR hijack*)',
}

EXCLUSION = (
    'NOT (movie* OR film* OR museum* OR anniversar* OR obituar* OR memorial* OR arts OR '
    'book OR books OR memoir* OR "price war" OR game OR story OR history OR veteran* OR '
    'tribute* OR sport OR music OR racing OR cancer OR "real estate" OR mafia OR trial OR tax)'
)


def _or_cats(ks):
    return " OR ".join(f"({CATEGORIES[k]})" for k in ks)


SERIES_KEYWORDS = {
    "GPR": _or_cats(range(1, 9)),
    "GPT": _or_cats(range(1, 6)),
    "GPA": _or_cats(range(6, 9)),
}
if CATEGORY_COUNTS:
    for _k in range(1, 9):
        SERIES_KEYWORDS[f"CAT{_k}"] = f"({CATEGORIES[_k]})"


def numerator_query(keywords):
    return f"({DTYPE_PREFIX} AND ({keywords}) {EXCLUSION})"


DENOM_QUERY = f'({DTYPE_PREFIX} AND ("THE" AND "BE" AND "TO" AND "OF" AND "AND" AND "AT" AND "IN"))'

# ---------------------------------------------------------------------------
# Corpus definitions
# ---------------------------------------------------------------------------


def _last_complete_month():
    t = date.today().replace(day=1) - timedelta(days=1)
    return (t.year, t.month)


def _flat(d, names=None):
    return [i for n, ids in d.items() if names is None or n in names for i in ids]


def paper_clauses(spec):
    """Query clauses selecting one paper: PUBID(...) per ID plus any title clauses."""
    return [f"PUBID({i})" for i in spec.get("pubids", [])] + list(spec.get("clauses", []))


def _all_clauses(papers, names=None):
    return [c for n, spec in papers.items() if names is None or n in names for c in paper_clauses(spec)]


def build_corpora():
    hist = {
        "key": "hist",
        "label": "Historical: NYT, Washington Post, Chicago Tribune",
        "products": list(HIST_PRODUCTS),
        "clauses": None,           # dedicated databases, no restriction needed
        "pubtitles": None,
        "start": (1900, 1),
        "end": (2022, 12),
        "ref": (1900, 2019),
        "bench": {"GPR": "GPRH", "GPT": "GPRHT", "GPA": "GPRHA"},
        "bench_N": "N3H", "bench_share": "SHARE_GPRH",
        "bench_cat_prefix": "SHAREH_CAT_",
        "probe_month": (1914, 8),  # published: N3H 7,542; GPR articles ~1,285
    }
    if HIST_EXTEND_WITH_CURRENT:
        three = ["New York Times", "Washington Post", "Chicago Tribune"]
        hist["products"] = list(HIST_PRODUCTS) + [p for p in RECENT_PRODUCTS if p not in HIST_PRODUCTS]
        hist["clauses"] = [f"PUBID({i})" for i in _flat(HIST_ARCHIVE_PUBIDS) + _flat(RECENT_PUBIDS, three)]
        hist["end"] = _last_complete_month()
        hist["label"] += " (archives + current editions)"

    recent = {
        "key": "recent",
        "label": ("Recent: C&I 10-paper panel (incl. Financial Times)" if INCLUDE_FT
                  else "Recent: C&I 10-paper panel minus FT (9 papers)"),
        "products": list(RECENT_PRODUCTS) + ([p for p in FT_PRODUCTS if p not in RECENT_PRODUCTS]
                                             if INCLUDE_FT else []),
        "clauses": _all_clauses(RECENT_PAPERS) if PUB_FILTER_MODE == "query" else None,
        "pubtitles": RECENT_PUBTITLES if PUB_FILTER_MODE == "pubtitle" else None,
        "start": (1985, 1),
        "end": _last_complete_month(),
        "ref": (1985, 2019),
        "bench": {"GPR": "GPR", "GPT": "GPRT", "GPA": "GPRA"},
        "bench_N": "N10", "bench_share": "SHARE_GPR",
        "bench_cat_prefix": None,
        "probe_month": (2001, 9),  # published: N10 33,230; SHARE_GPR 14.96%
        "papers": RECENT_PAPERS,
    }
    if PUB_FILTER_MODE == "pubtitle" and not RECENT_PUBTITLES:
        recent["_warn"] = ("PUB_FILTER_MODE='pubtitle' but RECENT_PUBTITLES is empty; "
                           "run `titles recent` and fill it in.")
    return {"hist": hist, "recent": recent}


# ---------------------------------------------------------------------------
# Low-level API (same as discover_facets.py)
# ---------------------------------------------------------------------------

def api_search(payload, max_retries=6):
    for attempt in range(1, max_retries + 1):
        try:
            r = requests.post(f"{BASE_URL}/api/cm/document/search", headers=headers,
                              json=payload, timeout=60)
            if r.status_code == 429:
                time.sleep(min(60, 2 ** attempt))
                continue
            r.raise_for_status()
            time.sleep(REQUEST_PAUSE)
            return r.json()
        except requests.exceptions.HTTPError as e:
            status = e.response.status_code if e.response is not None else None
            if status in (500, 502, 503, 504) and attempt < max_retries:
                time.sleep(min(30, 2 ** attempt))
                continue
            raise
        except requests.exceptions.RequestException:
            if attempt < max_retries:
                time.sleep(min(30, 2 ** attempt))
                continue
            raise
    raise RuntimeError("api_search: exhausted retries")


def _iso(d):
    return d.strftime("%Y-%m-%dT00:00:00.000Z")


def build_payload(corpus, query, start_iso, end_iso, facets=None, clauses=None):
    cl = clauses if clauses is not None else corpus.get("clauses")
    if cl:
        query = f"({query}) AND ({' OR '.join(cl)})"
    filters = [{"name": "sourcetype", "entries": [{"searchValue": "Newspapers"}]}]
    if corpus.get("pubtitles"):
        filters.append({"name": "pubTitle",
                        "entries": [{"searchValue": n} for n in corpus["pubtitles"]]})
    return {
        "workbenchId": WORKBENCH_ID,
        "search": {
            "query": query, "fulltext": True, "publications": [],
            "products": corpus["products"],
            "filters": filters,
            "startDate": start_iso, "endDate": end_iso,
        },
        "count": 1, "sortOrder": "pubDateDesc",
        "facets": [{"name": f} for f in (facets or [])],
    }


def get_count(corpus, query, start_iso, end_iso, clauses=None):
    return api_search(build_payload(corpus, query, start_iso, end_iso, clauses=clauses)).get("docsFound", 0)


# ---------------------------------------------------------------------------
# Month boundaries. discover_facets.py used endDate = 1st of next month; if the
# API treats endDate as inclusive that double-counts the 1st of each month.
# `detect_end_semantics` tests this once per corpus and caches the answer.
# ---------------------------------------------------------------------------

def month_bounds(year, month, end_inclusive):
    start = date(year, month, 1)
    nxt = date(year + (month == 12), month % 12 + 1, 1)
    end = nxt - timedelta(days=1) if end_inclusive else nxt
    return _iso(start), _iso(end)


def detect_end_semantics(corpus, verbose=True):
    os.makedirs(OUTPUT_ROOT, exist_ok=True)
    cache = os.path.join(OUTPUT_ROOT, f"_end_semantics_{corpus['key']}.json")
    if os.path.exists(cache):
        with open(cache) as f:
            return json.load(f)["end_inclusive"]
    y, m = corpus["probe_month"]
    day = _iso(date(y, m, 1))
    same_day = get_count(corpus, DENOM_QUERY, day, day)
    end_inclusive = same_day > 0
    # additivity check with the chosen convention
    a = get_count(corpus, DENOM_QUERY, *month_bounds(y, m, end_inclusive))
    y2, m2 = (y + (m == 12), m % 12 + 1)
    b = get_count(corpus, DENOM_QUERY, *month_bounds(y2, m2, end_inclusive))
    span = get_count(corpus, DENOM_QUERY, month_bounds(y, m, end_inclusive)[0],
                     month_bounds(y2, m2, end_inclusive)[1])
    if verbose:
        print(f"  [{corpus['key']}] start==end count on {day[:10]}: {same_day:,} -> "
              f"endDate treated as {'INCLUSIVE' if end_inclusive else 'EXCLUSIVE'}")
        print(f"  [{corpus['key']}] additivity: {a:,} + {b:,} = {a + b:,} vs two-month span {span:,} "
              f"{'OK' if a + b == span else '*** MISMATCH - check month boundaries ***'}")
    with open(cache, "w") as f:
        json.dump({"end_inclusive": end_inclusive, "same_day": same_day,
                   "a": a, "b": b, "span": span}, f)
    return end_inclusive


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def month_list(start, end):
    (y, m), out = start, []
    while (y, m) <= end:
        out.append((y, m))
        y, m = (y + (m == 12), m % 12 + 1)
    return out


def corpus_dir(corpus):
    d = os.path.join(OUTPUT_ROOT, corpus["key"])
    os.makedirs(d, exist_ok=True)
    return d


def collect(corpus):
    out_dir = corpus_dir(corpus)
    cp_path = os.path.join(out_dir, "checkpoint.json")
    cp = {}
    if os.path.exists(cp_path):
        with open(cp_path) as f:
            cp = json.load(f)
    lock = threading.Lock()
    end_inclusive = detect_end_semantics(corpus)
    needed = ["Denominator"] + list(SERIES_KEYWORDS)
    if "Financial Times" in corpus.get("papers", {}) and PUB_FILTER_MODE == "query":
        ft_counts = ft_coverage(corpus, end_inclusive, verbose=False)
        if not any(ft_counts.values()):
            msg = ("  *** Financial Times returned 0 documents in every test year: this run "
                   "is effectively the 9-paper index. Run `ft recent` to find FT. ***")
            if REQUIRE_FT:
                sys.exit(msg + "\n  (REQUIRE_FT = True, stopping.)")
            print(msg)
        else:
            print(f"  Financial Times found (denominator articles in test months: {ft_counts})")

    months = month_list(corpus["start"], corpus["end"])
    todo = [ym for ym in months
            if not all(k in cp.get(f"{ym[0]}-{ym[1]:02d}", {}) for k in needed)]
    print(f"\n[{corpus['key']}] {corpus['label']}")
    print(f"  {len(months)} months {months[0]}..{months[-1]}, {len(todo)} to fetch "
          f"({len(todo) * len(needed)} API calls)")

    def work(ym):
        key = f"{ym[0]}-{ym[1]:02d}"
        s, e = month_bounds(ym[0], ym[1], end_inclusive)
        row = dict(cp.get(key, {}))
        if "Denominator" not in row:
            row["Denominator"] = get_count(corpus, DENOM_QUERY, s, e)
        for name, kw in SERIES_KEYWORDS.items():
            if name not in row:
                row[name] = get_count(corpus, numerator_query(kw), s, e) if row["Denominator"] else 0
        return key, row

    done = 0
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futs = {ex.submit(work, ym): ym for ym in todo}
        for fut in as_completed(futs):
            try:
                key, row = fut.result()
            except Exception as e:
                print(f"  [ERROR] {futs[fut]}: {e}  (re-run to resume)")
                continue
            with lock:
                cp[key] = row
                done += 1
                if done % 25 == 0 or done == len(todo):
                    with open(cp_path, "w") as f:
                        json.dump(cp, f)
                    el = time.time() - t0
                    print(f"  [{corpus['key']}] {done}/{len(todo)} months "
                          f"({el / 60:.1f} min, ~{el / done * (len(todo) - done) / 60:.0f} min left)",
                          flush=True)
    with open(cp_path, "w") as f:
        json.dump(cp, f)
    missing = [ym for ym in months if f"{ym[0]}-{ym[1]:02d}" not in cp]
    if missing:
        print(f"  WARNING: {len(missing)} months still missing; re-run `run {corpus['key']}`.")
    return cp


# ---------------------------------------------------------------------------
# Index construction
# ---------------------------------------------------------------------------

def build_indices(corpus, cp):
    import pandas as pd
    rows = []
    for key, r in cp.items():
        y, m = map(int, key.split("-"))
        if not (corpus["start"] <= (y, m) <= corpus["end"]):
            continue
        rows.append({"month": pd.Timestamp(y, m, 1), **r})
    df = pd.DataFrame(rows).set_index("month").sort_index()
    den = df["Denominator"].where(df["Denominator"] > 0)
    r0, r1 = corpus["ref"]
    for name in SERIES_KEYWORDS:
        share = df[name] / den * 100           # percent, like C&I's SHARE_* columns
        df[f"SHARE_{name}"] = share
        if name in ("GPR", "GPT", "GPA"):
            df[f"{name}_index"] = share / share[str(r0):str(r1)].mean() * 100

    out_dir = corpus_dir(corpus)
    df.to_csv(os.path.join(out_dir, "monthly_counts.csv"), index_label="month")
    for name in ("GPR", "GPT", "GPA"):
        with open(os.path.join(out_dir, f"{name}_index.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Year", "Month", "Numerator", "Denominator", "Raw_Ratio",
                        f"{name}_Index_Normalized"])
            for t, r in df.iterrows():
                ok = r["Denominator"] > 0
                w.writerow([t.year, t.month, int(r[name]), int(r["Denominator"]),
                            r[name] / r["Denominator"] if ok else "",
                            r[f"{name}_index"] if ok else ""])
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.figure(figsize=(14, 6))
        for name, lw in (("GPR", 1.6), ("GPT", 1.2), ("GPA", 1.2)):
            plt.plot(df.index, df[f"{name}_index"], lw=lw, label=name, alpha=1 if name == "GPR" else 0.85)
        plt.title(f"{corpus['label']} (mean = 100 over {r0}-{r1})", fontsize=14)
        plt.legend(fontsize=12)
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, "indices.png"), dpi=150)
        plt.close()
    except ImportError:
        pass
    print(f"  [{corpus['key']}] wrote {out_dir}/GPR_index.csv, GPT_index.csv, GPA_index.csv")
    return df


# ---------------------------------------------------------------------------
# Validation against C&I's published file
# ---------------------------------------------------------------------------

def load_benchmark(path):
    import pandas as pd
    b = pd.read_csv(path) if path.lower().endswith(".csv") else pd.read_excel(path)
    # C&I's sheet starts Jan 1900; Excel's 1900 leap-year bug garbles the first
    # two date cells, so index by position and sanity-check the last row.
    idx = pd.date_range("1900-01-01", periods=len(b), freq="MS")
    try:
        last = pd.to_datetime(b["month"].iloc[-1])
        if (last.year, last.month) != (idx[-1].year, idx[-1].month):
            idx = pd.to_datetime(b["month"]).dt.to_period("M").dt.to_timestamp()
    except Exception:
        pass
    b.index = idx
    return b


def validate(corpus, df, bench_path):
    import pandas as pd
    bench = load_benchmark(bench_path)
    vdir = os.path.join(corpus_dir(corpus), "validation")
    os.makedirs(vdir, exist_ok=True)
    lines = [f"Validation: {corpus['label']}", f"Benchmark: {bench_path}", ""]

    comp = pd.DataFrame(index=df.index)
    for name, col in corpus["bench"].items():
        comp[f"{name}_ours"] = df[f"{name}_index"]
        comp[f"{name}_ci"] = bench[col].reindex(df.index)
    comp["Denominator_ours"] = df["Denominator"]
    comp[f"{corpus['bench_N']}_ci"] = bench[corpus["bench_N"]].reindex(df.index)
    comp["SHARE_GPR_ours"] = df["SHARE_GPR"]
    comp[f"{corpus['bench_share']}_ci"] = bench[corpus["bench_share"]].reindex(df.index)
    comp.to_csv(os.path.join(vdir, "comparison.csv"), index_label="month")

    lines.append("Correlation with published series (monthly):")
    res = {}
    for name, col in corpus["bench"].items():
        m = comp[[f"{name}_ours", f"{name}_ci"]].dropna()
        r = m.iloc[:, 0].corr(m.iloc[:, 1])
        res[name] = (r, len(m))
        lines.append(f"  {name} vs {col:6s}  r = {r:.4f}   (n = {len(m)}, "
                     f"{m.index.min():%Y-%m}..{m.index.max():%Y-%m})")

    # sub-periods
    lines += ["", "Sub-period correlations (GPR / GPT / GPA):"]
    edges = ([1900, 1920, 1946, 1970, 1985, 2000, 2010, 2020, 2100] if corpus["key"] == "hist"
             else [1985, 1987, 2000, 2010, 2020, 2100])
    for a, b_ in zip(edges[:-1], edges[1:]):
        sub = comp[str(a):str(b_ - 1)]
        rs = []
        for name in ("GPR", "GPT", "GPA"):
            s = sub[[f"{name}_ours", f"{name}_ci"]].dropna()
            rs.append(f"{s.iloc[:, 0].corr(s.iloc[:, 1]):.3f}" if len(s) > 2 else "  n/a")
        lab = f"{a}-{b_ - 1}" if b_ < 2100 else f"{a}-latest"
        lines.append(f"  {lab:12s} {'  '.join(rs)}   (n = {len(sub.dropna(subset=['GPR_ci']))})")

    # GPT-GPA correlation (C&I report 0.59 full sample, 0.45 from 1985)
    g = comp[["GPT_ours", "GPA_ours", "GPT_ci", "GPA_ci"]].dropna()
    lines += ["", f"GPT-GPA correlation: ours {g.GPT_ours.corr(g.GPA_ours):.3f}, "
                  f"C&I same window {g.GPT_ci.corr(g.GPA_ci):.3f}"]

    # level checks: shares and article counts
    s = comp[["SHARE_GPR_ours", f"{corpus['bench_share']}_ci"]].dropna()
    lines.append(f"Mean GPR share: ours {s.iloc[:, 0].mean():.2f}%  vs  C&I {s.iloc[:, 1].mean():.2f}%")
    py, pm = corpus["probe_month"]
    t = pd.Timestamp(py, pm, 1)
    if t in comp.index and pd.notna(comp.loc[t, f"{corpus['bench_N']}_ci"]):
        N = comp.loc[t, f"{corpus['bench_N']}_ci"]
        sh = comp.loc[t, f"{corpus['bench_share']}_ci"]
        lines.append(f"Single month {t:%Y-%m}: denominator ours {int(df.loc[t, 'Denominator']):,} vs "
                     f"C&I {int(N):,}; GPR articles ours {int(df.loc[t, 'GPR']):,} vs C&I "
                     f"~{round(N * sh / 100):,}")

    # optional category shares (historical benchmark has SHAREH_CAT_1..8)
    if CATEGORY_COUNTS and corpus.get("bench_cat_prefix"):
        lines += ["", "Category shares, correlation with C&I:"]
        for k in range(1, 9):
            col = f"{corpus['bench_cat_prefix']}{k}"
            if col in bench:
                x = pd.concat([df[f"SHARE_CAT{k}"], bench[col]], axis=1).dropna()
                lines.append(f"  cat {k}: r = {x.iloc[:, 0].corr(x.iloc[:, 1]):.3f}  "
                             f"mean share ours {x.iloc[:, 0].mean():.3f}% vs C&I {x.iloc[:, 1].mean():.3f}%")

    report = "\n".join(lines)
    print("\n" + report)
    with open(os.path.join(vdir, "report.txt"), "w") as f:
        f.write(report + "\n")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(3, 1, figsize=(14, 11), sharex=True)
        for ax, (name, col) in zip(axes, corpus["bench"].items()):
            m = comp[[f"{name}_ours", f"{name}_ci"]].dropna()
            ax.plot(m.index, m.iloc[:, 1], lw=0.9, label=f"C&I {col}", color="#999999")
            ax.plot(m.index, m.iloc[:, 0], lw=0.9, label="Replication", color="#990011", alpha=0.85)
            ax.set_title(f"{name}  (r = {res[name][0]:.3f}, n = {res[name][1]})", fontsize=13)
            ax.legend(loc="upper right")
            ax.grid(alpha=0.3)
        fig.suptitle(corpus["label"], fontsize=14)
        plt.tight_layout()
        plt.savefig(os.path.join(vdir, "overlay.png"), dpi=150)
        plt.close()
    except ImportError:
        pass
    return res


# ---------------------------------------------------------------------------
# Probe / titles helpers (run before a full build)
# ---------------------------------------------------------------------------

def probe(corpus):
    print(f"\n[probe] {corpus['label']}")
    print(f"  products: {[p['moniker'] for p in corpus['products']]}")
    if corpus.get("_warn"):
        print("  WARNING:", corpus["_warn"])
    end_inclusive = detect_end_semantics(corpus)
    y, m = corpus["probe_month"]
    s, e = month_bounds(y, m, end_inclusive)

    # Historical: test each HNP product on its own, so a bad moniker shows up as 0.
    if corpus["key"] == "hist" and len(corpus["products"]) > 1:
        print(f"  per-product denominators for {y}-{m:02d} (each should be > 0):")
        for prod in corpus["products"]:
            solo = dict(corpus)
            solo["products"] = [prod]
            c = get_count(solo, DENOM_QUERY, s, e)
            flag = "   <-- 0: WRONG MONIKER, open DevTools in the workbench and copy the real one" if c == 0 else ""
            print(f"    {prod['moniker']:30s} {c:>8,}{flag}")

    den = get_count(corpus, DENOM_QUERY, s, e)
    nums = {n: get_count(corpus, numerator_query(SERIES_KEYWORDS[n]), s, e) for n in ("GPR", "GPT", "GPA")}
    print(f"  {y}-{m:02d}: denominator {den:,}; GPR {nums['GPR']:,}  GPT {nums['GPT']:,}  GPA {nums['GPA']:,}"
          f"  -> GPR share {nums['GPR'] / den * 100 if den else float('nan'):.2f}%")
    if corpus["key"] == "hist":
        print("  expected (C&I, Aug 1914): denominator 7,542; GPR ~1,285 (share 17.04%). "
              "Build A got 7,516 / 1,274.")
    else:
        print("  C&I (Sep 2001, 10 papers): N10 33,230; SHARE_GPR 14.96%. Your counts will be larger "
              "(ProQuest indexes more items per paper), but the share should be in the same range.")
        print(f"  restriction clause: {' OR '.join(corpus.get('clauses') or [])[:200]}")
    if den == 0:
        print("  *** 0 documents: product moniker(s) or publication filter are wrong. "
              "See the PRODUCTS comment at the top of the script. ***")
    if nums["GPT"] + nums["GPA"] < nums["GPR"] or max(nums["GPT"], nums["GPA"]) > nums["GPR"]:
        print("  *** nesting check failed (expect GPT, GPA <= GPR <= GPT + GPA) ***")

    if corpus.get("papers") and PUB_FILTER_MODE == "query":
        print("  per-paper denominators (each should be > 0):")
        for paper, spec in corpus["papers"].items():
            cl = paper_clauses(spec)
            c = get_count(corpus, DENOM_QUERY, s, e, clauses=cl) if cl else 0
            ident = ",".join(spec.get("pubids") or []) or "title match"
            print(f"    {paper:22s} {ident:18s} {c:>8,}{'   <-- 0: check IDs/products' if c == 0 else ''}")
        if "Financial Times" in corpus["papers"]:
            ft_coverage(corpus, end_inclusive)


FT_TEST_YEARS = (1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2025)


def ft_coverage(corpus, end_inclusive, verbose=True):
    """Denominator-query count for the FT alone in one month (June) of several years."""
    cl = paper_clauses(corpus["papers"]["Financial Times"])
    out = {}
    for y in FT_TEST_YEARS:
        if (y, 6) > corpus["end"]:
            continue
        out[y] = get_count(corpus, DENOM_QUERY, *month_bounds(y, 6, end_inclusive), clauses=cl)
    if verbose:
        print("  Financial Times coverage (articles in June of each year):")
        print("    " + "  ".join(f"{y}: {c:,}" for y, c in out.items()))
        if not any(out.values()):
            print("    *** no FT documents: run `ft recent`, then set FT_PUBIDS / FT_TITLE_CLAUSES / "
                  "FT_PRODUCTS ***")
        elif any(c == 0 for c in out.values()):
            print("    FT coverage is partial; the same-corpus denominator absorbs this, as in C&I.")
    return out


def find_ft(corpus):
    """List every publication title ProQuest matches for the FT title clauses, by year,
    with no other paper restriction, so you can pick exact titles / spot false matches."""
    print(f"\n[ft] products searched: {[p['moniker'] for p in corpus['products']]}")
    print(f"     clauses: {FT_TITLE_CLAUSES}   ids: {FT_PUBIDS or '(none)'}")
    cl = paper_clauses(corpus["papers"]["Financial Times"]) if "Financial Times" in corpus.get("papers", {}) \
        else FT_TITLE_CLAUSES
    for y in FT_TEST_YEARS:
        s, e = _iso(date(y, 1, 1)), _iso(date(y, 12, 31))
        data = api_search(build_payload(corpus, DENOM_QUERY, s, e, facets=["pubTitle"], clauses=cl))
        block = next((b for b in data.get("facets", []) if b.get("name") == "pubTitle"), None)
        print(f"  {y}: {data.get('docsFound', 0):,} docs")
        for en in (block or {}).get("entries", [])[:10]:
            print(f"      {en.get('count', 0):>8,}  {en.get('displayValue', en.get('searchValue'))!r}")
    print("  If titles other than the FT appear, replace the clause with PUB.EXACT(\"<exact title>\").")
    print("  If nothing appears, FT is not in these products: add its database to FT_PRODUCTS.")


def search_pub_titles(corpus, pattern):
    """Scan the corpus for every publication whose title matches PUB(pattern), by year.
    Pass a short word: `financial`, `FT`, `times`, `tribune`. The pubTitle facet lists what's
    there under that pattern, so if FT exists under a different name (e.g. "Financial Times
    (London Edition)", "FT.com", "FT Weekend") it will show up here."""
    print(f"\n[search] pattern PUB({pattern}) across products {[p['moniker'] for p in corpus['products']]}")
    all_titles = {}
    for y in FT_TEST_YEARS:
        if (y, 1) > corpus["end"]:
            continue
        s_iso, e_iso = _iso(date(y, 1, 1)), _iso(date(y, 12, 31))
        # Minimal corpus for the search: no clause restriction, just DTYPE + pattern
        payload = {
            "workbenchId": WORKBENCH_ID,
            "search": {
                "query": f"{DTYPE_PREFIX} AND PUB({pattern})",
                "fulltext": True, "publications": [],
                "products": corpus["products"],
                "filters": [{"name": "sourcetype", "entries": [{"searchValue": "Newspapers"}]}],
                "startDate": s_iso, "endDate": e_iso,
            },
            "count": 1, "sortOrder": "pubDateDesc",
            "facets": [{"name": "pubTitle"}],
        }
        data = api_search(payload)
        print(f"  {y}: {data.get('docsFound', 0):,} docs matching PUB({pattern})")
        block = next((b for b in data.get("facets", []) if b.get("name") == "pubTitle"), None)
        for en in (block or {}).get("entries", [])[:15]:
            t = en.get("displayValue", en.get("searchValue", ""))
            c = en.get("count", 0)
            all_titles[t] = all_titles.get(t, 0) + c
            print(f"      {c:>8,}  {t!r}")
    if all_titles:
        print(f"\n  Distinct titles found (summed across years, top 20):")
        for t, c in sorted(all_titles.items(), key=lambda x: -x[1])[:20]:
            print(f"    {c:>10,}  {t!r}")
        print(f"\n  To use one of these, set:")
        print(f"    FT_TITLE_CLAUSES = ['PUB.EXACT(\"<pick a title above>\")']")
    else:
        print(f"\n  No publication titles matched PUB({pattern}) in any product in these products.")
        print(f"  Try the corpus-wide fallback: a facet scan with no PUB restriction, over 1 month,")
        print(f"  so the top-20 publications in the whole corpus show up — then look for FT.")
        print(f"  Or add another ProQuest product database to FT_PRODUCTS.")



def list_titles(corpus, year=None):
    y = year or corpus["probe_month"][0]
    s, e = _iso(date(y, 1, 1)), _iso(date(y, 12, 31))
    data = api_search(build_payload(corpus, DENOM_QUERY, s, e, facets=["pubTitle"]))
    block = next((b for b in data.get("facets", []) if b.get("name") == "pubTitle"), None)
    print(f"\n[titles] {corpus['label']}, {y}: docsFound {data.get('docsFound', 0):,}")
    for en in (block or {}).get("entries", []):
        print(f"  {en.get('count', 0):>9,}  {en.get('displayValue', en.get('searchValue'))!r}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["probe", "titles", "ft", "search", "run", "validate"])
    ap.add_argument("corpora", nargs="+", choices=["hist", "recent"])
    ap.add_argument("--bench", help="C&I data_gpr_export.xls (or .csv) for validation")
    ap.add_argument("--year", type=int, help="year for `titles`")
    ap.add_argument("--pattern", help="PUB(...) pattern for `search` (e.g. financial, FT, times)")
    args = ap.parse_args(argv)

    if TOKEN == "<TOKEN>" and args.command != "validate":
        sys.exit("Set TDM_TOKEN first:  export TDM_TOKEN='<bearer token>'")
    corpora = build_corpora()

    for key in args.corpora:
        c = corpora[key]
        if args.command == "probe":
            probe(c)
        elif args.command == "titles":
            list_titles(c, args.year)
        elif args.command == "ft":
            if key != "recent":
                print("`ft` applies to the recent corpus only.")
                continue
            find_ft(c)
        elif args.command == "search":
            search_pub_titles(c, args.pattern or "financial")
        elif args.command == "run":
            df = build_indices(c, collect(c))
            if args.bench:
                validate(c, df, args.bench)
        elif args.command == "validate":
            cp_path = os.path.join(OUTPUT_ROOT, key, "checkpoint.json")
            if not os.path.exists(cp_path):
                sys.exit(f"No data for {key}; run `run {key}` first.")
            with open(cp_path) as f:
                df = build_indices(c, json.load(f))
            if not args.bench:
                sys.exit("validate needs --bench")
            validate(c, df, args.bench)


if __name__ == "__main__":
    main()