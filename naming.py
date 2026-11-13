"""
naming.py  -  what the output tree calls things.

Separate from multimodel_etccdi.SCENARIO_LABEL on purpose. That table records
what the Reflective bucket calls each input directory, and the bucket is not
consistent: the SSP2-4.5 reference is 'SSP2-4.5' for CESM, 'SSP245' for UKESM
and 'baseline' for MIROC. Using it for the output made our own product inherit
someone else's inconsistency, which is what Alistair saw. The input table stays
exactly as it is, because it has to match the bucket to read it; this one says
what we write.

The scheme:

    <root>/<MODEL>/<SCENARIO>/<member>/
        <INDEX>_<MODEL>_<SCENARIO>_<member>_<first>-<last>.nc

for example

    MIROC-ES2H/SSP2-4.5/r01/TN10p_MIROC-ES2H_SSP2-4.5_r01_2020-2084.nc
    CESM2-WACCM/SSP2-4.5/r6/TX90p_CESM2-WACCM_SSP2-4.5_r6_2015-2069.nc
    UKESM1-1/G6-1.5K-SAI/r2/WSDI_UKESM1-1_G6-1.5K-SAI_r2_2035-2084.nc
    UKESM1-1/SSP2-4.5/r1/TN10p_UKESM1-1_SSP2-4.5_r1_2015-2100.nc
    UKESM1-0-LL/SSP2-4.5/r1/TN10p_UKESM1-0-LL_SSP2-4.5_r1_2015-2084.nc

Three choices in that, each of which can be changed by editing one table here
rather than hunting through the pipeline:

  the index comes first, so that a directory listing sorts by index. The
    directory path already carries model, scenario and member, so leading with
    those would sort every file identically and bury the one field that varies;

  the index keeps its ETCCDI capitalisation, TN10p rather than TN10P. The
    uppercase in the current tree comes from index_name.upper() in the writer
    and is not how the index is written anywhere in the literature;

  the year span IS in the name, following CMIP's habit, so a file carries its
    coverage without being opened. The cost is that extending a record renames
    its file, and a rerun that writes the longer span beside the shorter one
    leaves two copies of the same member. --overwrite on a rebuild, and
    verify_tree's extra-file report, are what catch that.

Member labels are not touched. They are r01 to r10 for MIROC and r1 to r3 or r6
to r10 for the others, which reflects how each archive numbers its own
realisations, and flattening them would break the correspondence with the
source data.
"""

# Output scenario directory names. Keys are the pipeline's internal scenario
# keys, not the bucket's directory names.
OUT_SCENARIO = {
    'SSP245': 'SSP2-4.5',
    'SAI':    'G6-1.5K-SAI',
    'HiLLA':  'G6-1.5K-HiLLA',
}

# Output model directory names.
#
# UKESM needs two entries, because the version differs by scenario rather than
# by model. The G6-1.5K-SAI and G6-1.5K-HiLLA runs are UKESM1-1 (Alistair, 18
# September 2026). The SSP2-4.5 reference exists in both versions: the one
# staged from the Reflective bucket is UKESM1-0-LL, confirmed from source_id
# and variant_label inside those files rather than from the directory name,
# which said UKESM1-1 and was wrong; the UKESM1-1 reference is the JASMIN
# macloud data.
#
# Until the second key existed, the tree carried one UKESM directory for all
# three scenarios, so whichever version it was named for mislabelled the other.
# baseline_version_note was carrying that distinction in an attribute, which no
# directory listing shows. Two model keys put it in the path instead, and both
# references can then sit in the tree at once, which is what lets the sensitivity
# of the comparison to the reference version be measured rather than argued.
OUT_MODEL = {
    'CESM':   'CESM2-WACCM',
    'UKESM':  'UKESM1-1',        # SAI, HiLLA, and the JASMIN SSP2-4.5
    'UKESM10': 'UKESM1-0-LL',    # the bucket SSP2-4.5 reference only
    'MIROC':  'MIROC-ES2H',
    'E3SM':   'E3SMv3',
}

# ETCCDI capitalisation, keyed by the lowercase form so a lookup works whatever
# case it is given.
INDEX_CANON = {
    'tn10p': 'TN10p', 'tn90p': 'TN90p', 'tx10p': 'TX10p', 'tx90p': 'TX90p',
    'r95p': 'R95p', 'r99p': 'R99p', 'wsdi': 'WSDI', 'csdi': 'CSDI',
    # The four fixed-threshold indices. RX1D and RX5D keep those exact names
    # rather than the ETCCDI spelling Rx1day and Rx5day, because that is what
    # the ReflectiveCloud/etccdi-indices zarr store calls them and the two index
    # sets are meant to merge without anyone renaming files. The ETCCDI
    # spellings map to the same output so a caller using either gets one file.
    'gsl': 'GSL', 'gdd': 'GDD',
    'rx1d': 'RX1D', 'rx5d': 'RX5D',
    'rx1day': 'RX1D', 'rx5day': 'RX5D',
    # The sixteen indices Cindy Wang and Alistair Duffey compute in
    # ReflectiveCloud/etccdi-indices. They are COMPUTED by the three
    # calculate_yearly_indices_*.ipynb notebooks, not by this pipeline; the
    # names live here because the output lands in the same member directories
    # and every module that walks the tree resolves names through this table.
    #
    # Each is spelled as their zarr store spells it, so the two sets merge
    # without anyone renaming a file. Two of those spellings are not the ETCCDI
    # ones and both are recorded rather than quietly corrected:
    #
    #   TN is their name for tropical nights, days with daily minimum
    #   temperature above 20 C. ETCCDI calls that index TR, and TN conventionally
    #   means daily minimum temperature itself, so this name will mislead a
    #   reader who knows the standard. 'tr' resolves here too, so code written
    #   either way finds the one file.
    #
    #   R10MM and R20MM are upper case in their store. This table previously
    #   mapped them to the ETCCDI spelling 'R10mm' and 'R20mm', which would have
    #   produced a second set of files beside theirs rather than the same set.
    'cdd': 'CDD', 'cwd': 'CWD', 'fd': 'FD', 'id': 'ID', 'su': 'SU',
    'tn': 'TN', 'tr': 'TN',
    'txx': 'TXX', 'txn': 'TXN', 'tnx': 'TNX', 'tnn': 'TNN',
    'prcptot': 'PRCPTOT', 'sdii': 'SDII',
    'r10mm': 'R10MM', 'r20mm': 'R20MM',
    'dtr': 'DTR',
}

# The span comes from the file's own time axis, so it states which years are
# present, not which are non-missing. UKESM SSP2-4.5 r2 spans 2015-2084 and has
# seven masked years inside that.
WITH_YEARS = True


# What each index is measured in, and the range a correct value lies in.
#
# This lives beside the canonical names because every module that reads the
# tree already imports this one, and because the alternative was a third copy:
# build_indices_derecho had its own INDEX_RANGE, verify_tree hardcoded 365 as
# "a full year", and the consistency notebook hardcoded the calendar length.
# All three were right while every index was a count of days and all three
# become wrong the moment one is not.
#
# The eight percentile indices and GSL count days. GDD is a sum of temperature
# differences in kelvin-days and reaches several thousand in the tropics. Rx1day
# and Rx5day are precipitation depths in millimetres. A check written against
# "days in the year" rejects the last three as broken when they are correct.
#
# The upper bounds are deliberately loose: they are there to catch a unit error
# or a corrupt write, which are wrong by orders of magnitude, not to police a
# plausible value. 15000 kelvin-days is about 41 K above the threshold every day
# of the year; 3000 mm in a day is above any measured total.
#
# RX5D is five times RX1D rather than a second guess. A five-day accumulation
# cannot exceed five times the largest daily total, so the bound follows from
# the daily one. The first value here was 6000, chosen by eye, and UKESM
# SSP2-4.5 reached 8527 mm, which rejected a correct file. That number is not
# physical at N96 resolution, the world record over five days being around
# 4800 mm, but it belongs to the same implausible precipitation tail already
# recorded for this model, where single days reach 1088 mm. A range check is
# not the tool for that: it catches a field that is wrong by a factor, and the
# tail is wrong by a fraction of a percent of the values.
# The strings here are the ones xclim stamps on the output, given this
# pipeline's inputs: temperature in kelvin and precipitation already converted
# to mm/day by derecho_io.read_member. They are not a tidied version of them.
# A table that disagrees with the file it describes is worse than no table, and
# this one did: GDD was recorded as 'K d' while every GDD file in the archive
# says 'd K'. Same quantity, and the file is what a user reads.
#
# The four annual extremes come back in the units of their input, which the
# xclim documentation states as "[same units as tasmax]". Our inputs are
# kelvin, so these files are in kelvin. The tables in the ReflectiveCloud
# notebooks say degrees Celsius for the same four indices, which is what those
# indices are conventionally reported in but not what the code produces, since
# it feeds kelvin straight through. Worth raising with them.
#
# Rx1day and Rx5day differ, 'mm d-1' against 'mm', and the reason is structural
# rather than cosmetic. max_1day_precipitation_amount resamples and takes a
# maximum, so it returns whatever it was given; max_n_day_precipitation_amount
# converts rate to amount first. Feeding the former kg m-2 s-1 returns a flux
# of about 1e-5 labelled kg m-2 s-1, which would sit inside any plausible
# range check and be wrong by nine orders of magnitude. This pipeline converts
# to mm/day before xclim sees it, as theirs does, so both are correct here.
INDEX_UNITS = {
    'TX90p': 'd', 'TX10p': 'd', 'TN90p': 'd', 'TN10p': 'd',
    'R95p': 'd', 'R99p': 'd', 'WSDI': 'd', 'CSDI': 'd',
    'GSL': 'd', 'GDD': 'd K', 'RX1D': 'mm d-1', 'RX5D': 'mm',
    # --- computed by the shared notebooks, not by this pipeline -------------
    # Read off xclim 0.54.0 by running each index, not copied from either
    # table. The four annual extremes come back in the units of their input,
    # which the xclim documentation states as "[same units as tasmax]"; the
    # notebooks feed kelvin, so these files are in kelvin. The tables inside the
    # ReflectiveCloud notebooks label the same four in Celsius, which is the
    # convention for reporting them but not what the code produces.
    'SU': 'd', 'ID': 'd', 'FD': 'd', 'TN': 'd',
    'CDD': 'd', 'CWD': 'd', 'R10MM': 'd', 'R20MM': 'd',
    'TXX': 'K', 'TXN': 'K', 'TNX': 'K', 'TNN': 'K',
    'PRCPTOT': 'mm', 'SDII': 'mm d-1',
}

INDEX_RANGE = {
    'TX90p': (0.0, 366.0), 'TX10p': (0.0, 366.0),
    'TN90p': (0.0, 366.0), 'TN10p': (0.0, 366.0),
    'R95p': (0.0, 366.0), 'R99p': (0.0, 366.0),
    'WSDI': (0.0, 366.0), 'CSDI': (0.0, 366.0),
    'GSL': (0.0, 366.0),
    'GDD': (0.0, 15000.0),
    'RX1D': (0.0, 3000.0), 'RX5D': (0.0, 5 * 3000.0),
    # --- computed by the shared notebooks ----------------------------------
    'SU': (0.0, 366.0), 'ID': (0.0, 366.0), 'FD': (0.0, 366.0),
    'TN': (0.0, 366.0), 'R10MM': (0.0, 366.0), 'R20MM': (0.0, 366.0),
    # CDD and CWD are longest spells, bounded by the year only because xclim
    # resamples before the run-length encoding. That is the default
    # (resample_before_rl=True) and what the notebooks use, so a multi-year
    # Saharan dry spell is cut at each new year rather than accumulating. If
    # that default ever changes these two can legitimately exceed a year and
    # this bound becomes wrong.
    'CDD': (0.0, 366.0), 'CWD': (0.0, 366.0),
    # Annual temperature extremes in kelvin. Impossibility bounds, not
    # plausibility ones: the coldest surface air temperature ever measured is
    # 184 K and the hottest 331 K, and UKESM's 1.5 m field in this archive
    # already reaches 367 K at a handful of desert cells, which is documented
    # and not corrected. What these catch is a field of zeros, or in Celsius,
    # or corrupt.
    'TXX': (150.0, 400.0), 'TXN': (150.0, 400.0),
    'TNX': (150.0, 400.0), 'TNN': (150.0, 400.0),
    # Annual total. The wettest place on Earth records about 26000 mm in a year
    # and a grid cell may exceed any observation, so the bound sits near four
    # times that. It catches a field left as a flux, not a wet cell.
    'PRCPTOT': (0.0, 100000.0),
    # Mean precipitation on wet days cannot exceed the largest daily total, so
    # it inherits the Rx1day bound rather than getting its own guess.
    'SDII': (0.0, 3000.0),
}

# Indices legitimately zero over large areas, where a high zero fraction is the
# climate rather than a truncated write. The spell indices are zero wherever no
# run of six days clears the threshold; GSL and GDD are zero wherever the
# temperature never reaches theirs, which is most of the high latitudes.
#
# The shared notebooks' counts are the same story with different geography. Ice
# days and frost days are zero throughout the tropics, summer days and tropical
# nights are zero throughout the high latitudes, and the two heavy-rain counts
# are zero over every desert. Each is most of a map, so a zero-fraction test
# would report the climate as a defect.
MOSTLY_ZERO_OK = {'WSDI', 'CSDI', 'GSL', 'GDD',
                  'SU', 'ID', 'FD', 'TN', 'R10MM', 'R20MM', 'CWD', 'CDD'}
# R99p was added here on 30 September 2026 and removed the same day. The
# reasoning for adding it was that the ETCCDI wet-day convention raises the
# threshold, so the annual count falls to about one percent of a cell's wet
# days, order one day a year, which would make an exactly-zero year ordinary.
# That was a prediction from the expected count and it was wrong. Measured over
# all 48 files of the pooled tree, R99p runs 20 to 32 percent zero across all
# three models and all three scenarios. R95p sits at 5 to 8 percent. WSDI and
# CSDI, at 34 to 77 percent, are why the set exists.
#
# Removing R99p from the set was still not the whole answer, because there are
# TWO guards and they trip at different places: build_indices_derecho's
# check_written above 0.5, verify_tree above 0.25. 20 to 32 percent clears the
# first and trips the second, so the first verify_tree run after the change
# reported 34 content problems, every one of them an R99p file behaving
# normally. A ceiling per index keeps both guards and drops the false alarms.


# Largest exactly-zero fraction that is climate rather than a truncated write,
# where the default for the calling guard is wrong. Set from measurement, not
# from expectation: R99p's observed maximum is 31.6 percent (MIROC-ES2H
# G6-1.5K-HiLLA r03, 30 September 2026), and 0.45 sits clear of that while
# staying well below what a partial write produces, since a truncated file is
# zero over a contiguous block rather than scattered.
ZERO_CEILING = {'R99p': 0.45}


def zero_ceiling(index_name, default):
    """
    The exactly-zero fraction above which this index is suspect.

    One place for a rule that two modules apply with different defaults.
    Membership of MOSTLY_ZERO_OK means no ceiling at all; otherwise a value in
    ZERO_CEILING overrides the caller's default, and anything else keeps it.
    """
    idx = canon_index(index_name)
    if idx in MOSTLY_ZERO_OK:
        return 1.0
    return ZERO_CEILING.get(idx, default)


# The two families this pipeline builds, in one place because four modules were
# each carrying their own copy of the list and they drifted: verify_tree and the
# consistency notebook disagreed about whether the four fixed-threshold indices
# were in the set for most of a day.
#
# PERCENTILE_INDICES need a threshold built from a reference run, so only they
# can be biased by the choice of reference. FIXED_INDICES do not.
#
# The other sixteen ETCCDI indices in ReflectiveCloud/etccdi-indices are NOT
# here. They are computed by the three calculate_yearly_indices_*.ipynb
# notebooks, which are the shared upstream code, and adding them to this module
# would mean two implementations of the same quantity to keep in step.
PERCENTILE_INDICES = ['TX90p', 'TX10p', 'TN90p', 'TN10p',
                      'R95p', 'R99p', 'WSDI', 'CSDI']
FIXED_INDICES = ['GSL', 'GDD', 'RX1D', 'RX5D']
ALL_INDICES = PERCENTILE_INDICES + FIXED_INDICES

# The sixteen from ReflectiveCloud/etccdi-indices, computed by the three
# calculate_yearly_indices_*.ipynb notebooks. Rx1day and Rx5day are in their
# set too and are already in FIXED_INDICES, so they are not repeated: this list
# is the fourteen the notebooks alone produce.
#
# TREE_INDICES is what a module that WALKS the tree should use, because both
# sets land in the same member directories. ALL_INDICES is what the build
# script uses, because it builds only its own twelve. Mixing the two is how a
# checker ends up reporting a perfectly good file as an extra.
SHARED_INDICES = ['SU', 'ID', 'TXX', 'TXN', 'FD', 'TN', 'TNX', 'TNN',
                  'PRCPTOT', 'SDII', 'CDD', 'CWD', 'R10MM', 'R20MM']
TREE_INDICES = ALL_INDICES + SHARED_INDICES


# The daily variable each index is built from, for all twenty-six.
#
# Here because a module that walks the tree needs to know which models can
# carry which index, and until now the only source of that was
# multimodel_etccdi.MODELS_FOR, which reads INDEX_REGISTRY. The registry holds
# the twelve this pipeline computes and by design does not hold the sixteen the
# shared notebooks compute, so MODELS_FOR raised KeyError on every one of them.
#
# That was not a loud failure. verify_tree caught the KeyError in a per-scenario
# except clause, printed "SSP245: 'SU'", skipped the scenario, and then reported
# '0 missing file(s), 0 content problem(s)' for an archive it had not looked at.
# A checker that passes because it checked nothing is worse than one that fails.
INDEX_VARIABLE = {
    'TX90p': 'tasmax', 'TX10p': 'tasmax', 'WSDI': 'tasmax',
    'TN90p': 'tasmin', 'TN10p': 'tasmin', 'CSDI': 'tasmin',
    'R95p': 'pr', 'R99p': 'pr', 'RX1D': 'pr', 'RX5D': 'pr',
    'GSL': 'tas', 'GDD': 'tas',
    'SU': 'tasmax', 'ID': 'tasmax', 'TXX': 'tasmax', 'TXN': 'tasmax',
    'FD': 'tasmin', 'TN': 'tasmin', 'TNX': 'tasmin', 'TNN': 'tasmin',
    'PRCPTOT': 'pr', 'SDII': 'pr', 'CDD': 'pr', 'CWD': 'pr',
    'R10MM': 'pr', 'R20MM': 'pr',
}

# Every model in this archive, and the ones whose daily maximum and minimum
# temperature are usable. E3SMv3 writes TREFHTMX == TREFHTMN == TREFHT at
# source, so an index built on either is meaningless for it.
MODELS_IN_TREE = ['CESM', 'UKESM', 'MIROC', 'E3SM']


def models_for(index_name):
    """
    Models that can carry this index, for any of the twenty-six.

    Unknown indices return every model rather than raising, so a name this table
    has not caught up with is reported as missing where it might exist instead
    of removing a whole scenario from a check.
    """
    var = INDEX_VARIABLE.get(canon_index(index_name))
    if var in ('tasmax', 'tasmin'):
        return [m for m in MODELS_IN_TREE if m != 'E3SM']
    return list(MODELS_IN_TREE)


# Days in a year, by calendar. Here rather than in baseline_check because
# verify_tree needs it too, and its absence there was a live defect: the test
# for an index pinned at the length of the year compared against a hardcoded
# 365, which no value on UKESM's 360-day calendar can reach, so the check that
# catches a corrupt threshold was inert for that model on every index.
DAYS_IN_YEAR = {
    '360_day': 360.0,
    'noleap': 365.0, '365_day': 365.0,
    'all_leap': 366.0, '366_day': 366.0,
    'standard': 365.2425, 'gregorian': 365.2425,
    'proleptic_gregorian': 365.2425, 'julian': 365.25,
}

# cftime class names, for the case where .dt.calendar is unavailable.
CFTIME_CALENDAR = {
    'Datetime360Day': '360_day', 'DatetimeNoLeap': 'noleap',
    'DatetimeAllLeap': 'all_leap', 'DatetimeJulian': 'julian',
    'DatetimeGregorian': 'standard',
    'DatetimeProlepticGregorian': 'proleptic_gregorian',
}


def calendar_of(da):
    """
    The calendar of a DataArray's time axis, or None.

    Three places carry it and any one of them can be missing depending on how
    the file was written and decoded, so all three are tried rather than
    assuming the first.
    """
    try:
        cal = da['time'].dt.calendar
        if cal:
            return str(cal)
    except (AttributeError, TypeError, ValueError, KeyError):
        pass
    for src in (getattr(da['time'], 'encoding', {}),
                getattr(da['time'], 'attrs', {})):
        if src.get('calendar'):
            return str(src['calendar'])
    try:
        return CFTIME_CALENDAR.get(type(da['time'].values[0]).__name__)
    except (IndexError, TypeError, KeyError):
        return None


def days_in_year(calendar, default=366.0):
    """
    Length of a year on this calendar.

    An unknown calendar returns the default of 366, the longest ordinary year,
    which makes a saturation test conservative rather than spuriously strict.
    """
    if calendar is None:
        return default
    return DAYS_IN_YEAR.get(str(calendar), default)


def max_days_in_year(calendar, default=366.0):
    """
    Days in the LONGEST year on this calendar.

    Distinct from days_in_year, and the distinction matters. DAYS_IN_YEAR holds
    365.2425 for a Gregorian calendar, which is the right number for "a tenth
    of a year" in the percentile check and the wrong one for "the length of the
    year" in a saturation check: growing season length reaches 366 in a leap
    year, which is correct and which a test written against 365.2425 reports as
    over-long. Every calendar with leap days answers 366 here; the fixed ones
    answer their own length.
    """
    if calendar is None:
        return default
    n = DAYS_IN_YEAR.get(str(calendar))
    if n is None:
        return default
    # A fractional length means the calendar has leap years, so the longest
    # year is the next whole number up.
    return float(int(n)) if float(n).is_integer() else float(int(n) + 1)


def days_in_calendar_year(year, calendar, default=None):
    """
    Days in ONE named year on this calendar, exactly.

    The third member of the family, and the one the other two cannot stand in
    for. days_in_year answers 365.2425 for a Gregorian calendar, the mean year,
    which is right for "a tenth of a year" and wrong for any particular one.
    max_days_in_year answers 366, the longest, which is right for a saturation
    bound and wrong as an expectation: a check that expects 366 of every
    Gregorian year reports all 37 common years in a 50-year record as a day
    short. That happened, in check_spans, and it buried three members that
    genuinely stop part-way through 2084 under 370 lines of false alarm.

    Returns None for a calendar it does not know rather than guessing, so a
    caller can say it does not know instead of reporting a difference.
    """
    cal = str(calendar) if calendar is not None else ''
    fixed = {'360_day': 360, 'noleap': 365, '365_day': 365,
             'all_leap': 366, '366_day': 366}
    if cal in fixed:
        return fixed[cal]
    y = int(year)
    if cal == 'julian':
        return 366 if y % 4 == 0 else 365
    if cal in ('proleptic_gregorian', 'standard', 'gregorian'):
        # The mixed Julian/Gregorian 'standard' calendar differs from the
        # proleptic one only before 1582, and nothing in this archive reaches
        # back that far, so one rule covers all three here.
        return 366 if (y % 4 == 0 and y % 100 != 0) or y % 400 == 0 else 365
    return default


def is_day_count(index_name):
    """
    True when the index is a number of days, so a bound of one year applies.

    Unknown indices return False, which disables the day-count checks rather
    than applying a bound that may be meaningless. A check that does not run is
    visible in the output; one that runs against the wrong scale reports a
    defect that is not there.
    """
    return INDEX_UNITS.get(canon_index(index_name)) == 'd'


def index_range(index_name):
    """
    (low, high) a correct value lies in, or (-inf, inf) if not recorded.
    """
    return INDEX_RANGE.get(canon_index(index_name), (float('-inf'),
                                                     float('inf')))


def out_scenario(scenario):
    if scenario not in OUT_SCENARIO:
        raise KeyError(f'no output name for scenario {scenario!r}; '
                       f'known: {sorted(OUT_SCENARIO)}')
    return OUT_SCENARIO[scenario]


def out_model(model):
    if model not in OUT_MODEL:
        raise KeyError(f'no output name for model {model!r}; '
                       f'known: {sorted(OUT_MODEL)}')
    return OUT_MODEL[model]


def canon_index(index_name):
    """
    ETCCDI capitalisation. Unknown indices pass through unchanged rather
    than being forced to a guess.
    """
    return INDEX_CANON.get(index_name.lower(), index_name)


def index_filename(index_name, model, scenario, member_label,
                   years=None, with_years=None):
    """
    The filename for one index of one member.

    years is (first, last) and is used only when years are being included.
    """
    use = WITH_YEARS if with_years is None else with_years
    parts = [canon_index(index_name), out_model(model), out_scenario(scenario),
             member_label]
    if use and years:
        parts.append(f'{years[0]}-{years[1]}')
    return '_'.join(parts) + '.nc'


def find_index_file(member_dir, index_name):
    """
    Path to one index inside a member directory, under either scheme.

    Readers go through this so the tree can be renamed without every script
    needing to change on the same day, and so a half-renamed tree, which is
    what an interrupted rename leaves, still reads.

    Returns (path, style) with style 'new' or 'old', or (None, None).

    Raises when more than one new-style file in the directory claims the same
    index. That is not a hypothetical: on 19 September 2026 a UKESM rebuild
    wrote TN10p_UKESM1-1_... beside a TN10p_UKESM1-0-LL_... that a pending
    rename had never removed, and this function returned the older one because
    '0' sorts before '1'. verify_tree reported the duplicates as extra files;
    analyze_tree, which comes through here, silently read the stale set and
    produced a table comparing indices built from two different percentile
    thresholds. A reader cannot be expected to notice that in a number.

    So an ambiguous directory is an error rather than a choice. Picking the
    newest by modification time would be worse, because it would look like it
    worked.
    """
    import os
    want = canon_index(index_name).lower()
    old_names = {f'{index_name.upper()}.nc', f'{canon_index(index_name)}.nc'}
    matches, fallback = [], None
    try:
        entries = os.listdir(member_dir)
    except OSError:
        return None, None
    for fn in sorted(entries):
        if not fn.endswith('.nc'):
            continue
        parsed = parse_index_filename(fn)
        if parsed and parsed['index'].lower() == want:
            matches.append(fn)
        elif fn in old_names:
            fallback = os.path.join(member_dir, fn)
    if len(matches) > 1:
        raise RuntimeError(
            f'{member_dir}: {len(matches)} files claim index {index_name!r}: '
            f'{matches}. Two builds have written into the same member '
            f'directory, or a rename was planned and never applied. Remove the '
            f'ones you do not want before reading anything from this tree; '
            f'whichever this returned would be a guess.')
    if matches:
        return os.path.join(member_dir, matches[0]), 'new'
    return (fallback, 'old') if fallback else (None, None)


def parse_index_filename(name):
    """
    Pull the fields back out of a name this module made.

    Returns a dict, or None when the name does not match the scheme, which is
    how a reader tells a new-style file from the old bare INDEX.nc.
    """
    if not name.endswith('.nc'):
        return None
    stem = name[:-3]
    bits = stem.split('_')
    if len(bits) < 4:
        return None
    out = {'index': bits[0], 'model': bits[1], 'scenario': bits[2],
           'member': bits[3], 'years': None}
    if len(bits) >= 5 and '-' in bits[4]:
        a, _, b = bits[4].partition('-')
        if a.isdigit() and b.isdigit():
            out['years'] = (int(a), int(b))
    return out