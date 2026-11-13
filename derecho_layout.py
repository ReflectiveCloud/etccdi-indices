"""
derecho_layout.py  -  where the daily fields live on GLADE.

Every path in this file is transcribed from the team spreadsheet
SAI_extremes_DataPaper_Variables.xlsx, sheet 'variables', as received on
14 September 2026. I have not listed any of these directories on Derecho, so
they are the team's record, not a verified listing. Check one file per
(model, scenario) before a production run; the spreadsheet is a planning
document and paths in it may predate a move.

Shape difference from the bucket. MODEL_BUCKET_LAYOUT in multimodel_etccdi.py
describes directories to be globbed. The GLADE tree is described per variable
with a filename template and an {ensemble_id} placeholder, and the template
differs by variable within a scenario, so the layout here is keyed
(model, scenario) -> {var: template}. The reader substitutes ensemble_id and
globs; templates that already end in a directory are globbed with '*.nc'.

Gaps are left as None and raise on access. Two of them matter:

  CESM SSP2-4.5 was missing from the spreadsheet and was supplied by Cindy on
  16 September 2026: /glade/campaign/cesm/collections/CESM2-WACCM-SSP245,
  members 006-010 only, because 1-5 carry a bug in daily minimum and maximum
  temperature. Confirmed by listing the collection, which quarantines the broken
  files in TREFHTMN.h1.bad and TREFHTMX.h1.bad.

  E3SMv3 has no rows at all in the spreadsheet. E3SM carries no temperature
  indices anyway, so this blocks only R95p, R99p and the fixed-threshold
  precipitation indices for that model. It also raises the question of whether
  E3SM is in the Derecho build at all, or whether it stays on the hub.
"""

import os

# --------------------------------------------------------------------------
# Member sets. These are the spreadsheet's counts, which differ from the
# Reflective bucket for MIROC.
# --------------------------------------------------------------------------
MEMBERS = {
    ('CESM',  'SAI'):    ('001', '002', '003'),         # count 3; ids to confirm
    ('CESM',  'HiLLA'):  ('001', '002', '003'),         # count 3; ids to confirm
    # Confirmed on disk 16 September 2026. Directory ids are zero-padded, and
    # only 006-010 are used: Cindy reported that members 1-5 have a bug in daily
    # minimum and maximum temperature, and the archive agrees with her. The
    # broken files were quarantined into TREFHTMN.h1.bad and TREFHTMX.h1.bad at
    # the top of the collection, and the ones listed there belong to members 001
    # and 002 onward.
    ('CESM',  'SSP245'): ('006', '007', '008', '009', '010'),

    # MIROC is 10 members on GLADE, not the 3 staged on the Reflective bucket.
    ('MIROC', 'SAI'):    tuple(f'{i:02d}' for i in range(1, 11)),
    ('MIROC', 'HiLLA'):  tuple(f'{i:02d}' for i in range(1, 11)),
    ('MIROC', 'SSP245'): tuple(f'{i:02d}' for i in range(1, 11)),

    # Listed on disk 15 September 2026 (01_followup_diagnostics section A).
    # They are plain 001/002/003, not Met Office u-XXXXX suite ids as the
    # filename template suggested.
    ('UKESM', 'SAI'):    ('001', '002', '003'),
    ('UKESM', 'HiLLA'):  ('001', '002', '003'),

    # UKESM SSP2-4.5 now has two entries, because two model versions of the
    # reference exist and the scenarios are UKESM1-1. 'UKESM' is the UKESM1-1
    # reference fetched from the JASMIN macloud share on 19 September 2026;
    # 'UKESM10' is the UKESM1-0-LL reference already staged on GLADE, kept so
    # the comparison's sensitivity to the reference version can be measured.
    #
    # The three JASMIN suites carry no variant_label, realization_index or
    # source_id: their only global attribute of any use is 'Data from Met
    # Office Unified Model'. So r1, r2 and r3 below are OT's provisional
    # numbering from the suite ids (19 September 2026), pending confirmation
    # from Alistair, and not a claim that these are r1i1p1f2 and so on. The
    # provenance is written into each output file by build_indices_derecho so
    # a reader is not left to infer it from a directory name.
    ('UKESM',   'SSP245'): ('u-dp902', 'u-dp903', 'u-dp904'),
    ('UKESM10', 'SSP245'): ('001', '002', '003'),
}

# Last year present in each member's files, where it is not the scenario norm.
# From 01_followup_diagnostics section D: MIROC HiLLA r04 and r06 are genuinely
# shorter runs, not truncated transfers. File size scales with record length
# (1.51 GB over 43 years and 1.58 GB over 45, against 1.75 GB over 50) and every
# file carries the same 25 August 2026 modification time, so a partial sync is
# ruled out. Both fall short of the 2065-2084 evaluation window.
SHORT_MEMBERS = {
    ('MIROC', 'HiLLA'): {},
}

# Members excluded from the percentile baseline, with the reason.
#
# UKESM 002 (variant r2i1p1f2) is missing seven whole years inside 2020-2039 on
# GLADE: 2020, 2021, 2022, 2025, 2026, 2030 and 2035, in three separate gaps,
# leaving 13 years of 20. Confirmed for both tasmin and pr, so it is the member
# rather than one variable. A threshold averaged over members of unequal length
# weights them unequally for no physical reason, and 13 years gives 65 samples
# per calendar day against 100.
#
# The Reflective bucket independently uses realizations 12 and 3 for UKESM
# SSP2-4.5 and not realization 2, which is at least consistent with r2 being
# known to be incomplete.
#
# This belongs to UKESM10, the GLADE UKESM1-0-LL reference, not to the JASMIN
# UKESM1-1 reference, whose three suites each carry a complete 2015-2100 record
# of 30,960 daily steps on a 360_day calendar with no gaps (probe_ukesm11.py,
# 19 September 2026). Leaving this keyed on 'UKESM' would have excluded a
# UKESM1-1 member for a defect in a different model version's file.
BASELINE_EXCLUDE = {
    ('UKESM10', 'SSP245'): {'002': 'missing 7 of 20 years in 2020-2039 '
                                   '(2020-2022, 2025, 2026, 2030, 2035)'},
}

# The reference ensemble is cut to the size of the SAI ensemble, agreed with
# Cindy Wang on 30 September 2026. Her wording: use the same number of SSP2-4.5
# members as the G6 experiments have, and for CESM2-WACCM she uses 006, 007 and
# 008, so those are the three named here rather than any three of the five.
#
# It binds on CESM2-WACCM alone. UKESM1-1 is already 3 against 3. MIROC-ES2H is
# 10 against 10 and stays at 10: its HiLLA ensemble IS ten members, of which r04
# and r06 stop in 2077 and 2079 (SHORT_MEMBERS above) and so drop out of the
# 2065-2084 evaluation window. Matching the baseline to the 8 that survive that
# window would give MIROC one threshold for HiLLA and another for SAI, and
# HiLLA minus SAI is the contrast most of the results rest on, so the two arms
# have to share a threshold.
#
# The cost of the cut was measured before it was made, on CESM over 30N to 60N:
# the 5-member and 3-member pooled thresholds differ by a median of -0.0002
# mm/day for R99p and -0.035 K for TX90p, with an rms of 6.45 mm/day and 0.30 K.
# No systematic shift, added scatter (plot_threshold_conventions.py, panel d).
BASELINE_MATCH_SAI = {
    ('CESM', 'SSP245'): ('006', '007', '008'),
}


# Directory id on GLADE -> member label used in the output tree, where the two
# differ. CESM SSP2-4.5 directories are zero-padded ('006') while the Reflective
# bucket writes 'r6', and the two trees should not disagree on what a member is
# called.
_MIROC_RNN = {f'{i:02d}': f'r{i:02d}' for i in range(1, 11)}

OUT_LABEL = {
    ('CESM',  'SSP245'): {'006': 'r6', '007': 'r7', '008': 'r8',
                          '009': 'r9', '010': 'r10'},
    ('CESM',  'SAI'):    {'001': 'r1', '002': 'r2', '003': 'r3'},
    ('CESM',  'HiLLA'):  {'001': 'r1', '002': 'r2', '003': 'r3'},
    
    # The bucket writes MIROC members as r01, r02, r03. GLADE directories are
    # bare '01' to '10', so without this the two trees disagree on the member
    # name for MIROC and a reader cannot tell 01 from r01 without opening a file.
    ('MIROC', 'SSP245'): _MIROC_RNN,
    ('MIROC', 'HiLLA'):  _MIROC_RNN,
    ('MIROC', 'SAI'):    _MIROC_RNN,
    # UKESM SAI maps through MEMBER_LABEL ('001' -> 'r1'), but HiLLA's
    # MEMBER_LABEL is keyed on the suite ids the bucket uses (r12i1p1f2 and so
    # on), which GLADE ids do not match, so HiLLA would fall through unmapped
    # and write '001' beside SAI's 'r1'. Mapped here so the tree is uniform.
    ('UKESM', 'HiLLA'):  {'001': 'r1', '002': 'r2', '003': 'r3'},
    ('UKESM', 'SAI'):    {'001': 'r1', '002': 'r2', '003': 'r3'},

    # Confirmed by Alistair Duffey on 25 September 2026 from the team
    # spreadsheet: u-dp902 is r2i1p1f2, u-dp903 is r3i1p1f2 and u-dp904 is
    # r12i1p1f2. Our earlier labels r1, r2, r3 came from the numerical order of
    # the suite ids, because the JASMIN files carry no variant_label, source_id
    # or realization_index at all, and that ordering was wrong for all three.
    # The written tree was corrected with src/relabel_members.py the same day.
    ('UKESM',   'SSP245'): {'u-dp902': 'r2', 'u-dp903': 'r3',
                            'u-dp904': 'r12'},
    ('UKESM10', 'SSP245'): {'001': 'r1', '002': 'r2', '003': 'r3'},
}


def out_label(model, scenario, member):
    """Member directory name in the output tree.

    Raises on a member with no mapping that still looks like a raw GLADE id,
    rather than writing a directory whose name means something different from
    its neighbours. A silently unmapped member is how a tree ends up with '001'
    and 'r1' side by side.
    """
    table = OUT_LABEL.get((model, scenario))
    if table is None:
        return member
    if member not in table:
        raise KeyError(
            f'no output label for {model}/{scenario}/{member}. Add it to '
            f'OUT_LABEL; leaving it unmapped would write a directory name '
            f'inconsistent with the rest of the tree.')
    return table[member]


def baseline_members(model, scenario='SSP245'):
    """
    Members the percentile threshold is pooled over, with the excluded ones
    dropped and named.

    Two filters, in order. BASELINE_EXCLUDE removes members with defective
    records. BASELINE_MATCH_SAI then cuts what is left to the size of the SAI
    ensemble, which is a convention rather than a data problem.
    """
    all_m = members(model, scenario)
    drop = BASELINE_EXCLUDE.get((model, scenario), {})
    keep = tuple(m for m in all_m if m not in drop)
    for m, why in drop.items():
        if m in all_m:
            print(f'  baseline excludes {model}/{scenario}/{m}: {why}')
    if not keep:
        raise RuntimeError(f'every baseline member excluded for {model}')

    match = BASELINE_MATCH_SAI.get((model, scenario))
    if match:
        missing = [m for m in match if m not in keep]
        if missing:
            # Naming a member that is not there means the archive moved under
            # the convention. Guessing a substitute would put an unrecorded
            # member set into every threshold, so stop instead.
            raise RuntimeError(
                f'BASELINE_MATCH_SAI names {missing} for {model}/{scenario}, '
                f'which is not among the available members {keep}. Fix the '
                f'convention or the layout; do not substitute.')
        print(f'  baseline cut to {len(match)} of {len(keep)} member(s) to '
              f'match the SAI ensemble size: {", ".join(match)}')
        keep = tuple(match)
    return keep

_MIROC = '/glade/derecho/scratch/cindywang625/other_models/MIROC'
_UKESM = '/glade/derecho/scratch/cindywang625/other_models/UKESM'

# UKESM1-1 SSP2-4.5, fetched from
# https://gws-access.jasmin.ac.uk/public/macloud/SSP245_ukesm1p1/ on
# 19 September 2026 by fetch_ukesm11_ssp245.py. This one IS a verified listing:
# every path below was read off the server and the temperature files were
# opened and checked (probe_ukesm11.py). Set REFLECTIVE_UKESM11_DIR to move it.
_UKESM11 = os.environ.get(
    'REFLECTIVE_UKESM11_DIR',
    os.path.join('/glade/work', os.environ.get('USER', 'nobody'),
                 'ukesm1p1_ssp245'))
# SSP2-4.5 lives in its own collection, not beside the scenarios, and its case
# name differs: BWSSP245cmip6 rather than BW...SSP245-G6-1p5K-*. The daily
# stream is day_1, the same as the scenarios. Member 006 also has a companion
# '006ext' directory holding its post-2069 extension; the baseline ends in 2034
# so the base member covers it and the extension is not read.
_CESM_SSP245 = ('/glade/campaign/cesm/collections/CESM2-WACCM-SSP245/'
                'b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.{ensemble_id}/'
                'atm/proc/tseries/day_1/'
                'b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.'
                '{ensemble_id}.cam.h1.{cam_var}.*.nc')

_CESM_SAI = ('/glade/campaign/cesm/collections/ARISE-SAI-1.5/'
             'b.e21.BW.f09_g17.SSP245-G6-1p5K-SAI.{ensemble_id}/atm/proc/'
             'tseries/day_1/b.e21.BW.f09_g17.SSP245-G6-1p5K-SAI.'
             '{ensemble_id}.cam.h1.{cam_var}.*.nc')
_CESM_HILLA = ('/glade/campaign/cgd/amp/walkerl/GeoMIP/HiLLA/'
               'b.e21.BW.f09_g17.SSP245-G6-1p5K-HiLLA.{ensemble_id}/atm/proc/'
               'tseries/day_1/b.e21.BW.f09_g17.SSP245-G6-1p5K-HiLLA.'
               '{ensemble_id}.cam.h1.{cam_var}.*.nc')

CAM_VAR = {'pr': 'PRECT', 'tas': 'TREFHT',
           'tasmax': 'TREFHTMX', 'tasmin': 'TREFHTMN'}

LAYOUT = {
    # ---- CESM2-WACCM -------------------------------------------------------
    ('CESM', 'SAI'):   {v: _CESM_SAI.replace('{cam_var}', c)
                        for v, c in CAM_VAR.items()},
    ('CESM', 'HiLLA'): {v: _CESM_HILLA.replace('{cam_var}', c)
                        for v, c in CAM_VAR.items()},
    ('CESM', 'SSP245'): {v: _CESM_SSP245.replace('{cam_var}', c)
                         for v, c in CAM_VAR.items()},

    # ---- MIROC-ES2H --------------------------------------------------------
    ('MIROC', 'SAI'): {
        'pr':     f'{_MIROC}/G6-1.5K-SAI/pr/pr_G6-1.5K-SAI_r{{ensemble_id}}.nc',
        'tas':    f'{_MIROC}/G6-1.5K-SAI/tas/tas_G6-1.5K-SAI_r{{ensemble_id}}.nc',
        'tasmax': f'{_MIROC}/G6-1.5K-SAI/tasmax/tasmax_G6-1.5K-SAI_r{{ensemble_id}}.nc',
        'tasmin': f'{_MIROC}/G6-1.5K-SAI/tasmin/tasmin_G6-1.5K-SAI_r{{ensemble_id}}.nc',
    },
    ('MIROC', 'HiLLA'): {
        'pr':     f'{_MIROC}/HiLLA-SAI/daily/pr/pr_G6-1.5K-HiLLA_r{{ensemble_id}}.nc',
        'tas':    f'{_MIROC}/HiLLA-SAI/daily/tas/tas_G6-1.5K-HiLLA_r{{ensemble_id}}.nc',
        'tasmax': f'{_MIROC}/HiLLA-SAI/daily/tasmax/tasmax_G6-1.5K-HiLLA_r{{ensemble_id}}.nc',
        'tasmin': f'{_MIROC}/HiLLA-SAI/daily/tasmin/tasmin_G6-1.5K-HiLLA_r{{ensemble_id}}.nc',
    },
    ('MIROC', 'SSP245'): {
        'pr':     f'{_MIROC}/SSP2-4.5/pr/pr_baseline_r{{ensemble_id}}.nc',
        'tas':    f'{_MIROC}/SSP2-4.5/tas/tas_baseline_r{{ensemble_id}}.nc',
        'tasmax': f'{_MIROC}/SSP2-4.5/tasmax/tasmax_baseline_r{{ensemble_id}}.nc',
        'tasmin': f'{_MIROC}/SSP2-4.5/tasmin/tasmin_baseline_r{{ensemble_id}}.nc',
    },

    # ---- UKESM1-1 ----------------------------------------------------------
    # G6-1.5K-SAI is the awkward one: tas, tasmax and tasmin all sit in the same
    # file, distinguished only by internal variable name. See VAR_IN_FILE below.
    # Confirmed on disk: this scenario has pr/ and tas/ directories only. There
    # are no tasmax/ or tasmin/ directories, because all three temperature
    # fields sit inside the tas file under different internal names. See
    # VAR_IN_FILE below, and verify that mapping before use.
    # Listed on disk 30 September 2026: daily/ holds pr, tas and rsdsdiff, and
    # each variable directory holds one directory per member, so the file sits
    # at daily/pr/001/pr_day_UKESM1-1-LL_u-001_201501-210012.nc. An earlier
    # version of this entry named the file directly under daily/pr/ and omitted
    # the member directory, so the glob matched nothing and the scenario read as
    # absent. Trailing slashes are load-bearing: the reader appends '*.nc' only
    # to a path that ends in one.
    ('UKESM', 'SAI'): {
        'pr':     f'{_UKESM}/G6-1.5K-SAI/daily/pr/{{ensemble_id}}/',
        'tas':    f'{_UKESM}/G6-1.5K-SAI/daily/tas/{{ensemble_id}}/',
        'tasmax': f'{_UKESM}/G6-1.5K-SAI/daily/tas/{{ensemble_id}}/',
        'tasmin': f'{_UKESM}/G6-1.5K-SAI/daily/tas/{{ensemble_id}}/',
    },
    ('UKESM', 'HiLLA'): {
        'pr':     f'{_UKESM}/HiLLA-SAI/daily/pr/{{ensemble_id}}/',
        'tas':    f'{_UKESM}/HiLLA-SAI/daily/tas/{{ensemble_id}}/',
        # Trailing slashes matter: the reader appends '*.nc' only to a path
        # that ends in one, so without them the glob returns the directory
        # itself and open_dataset fails with an engine error that looks like a
        # missing netCDF backend rather than a path bug.
        'tasmax': f'{_UKESM}/HiLLA-SAI/daily/tasmax/{{ensemble_id}}/',
        'tasmin': f'{_UKESM}/HiLLA-SAI/daily/tasmin/{{ensemble_id}}/',
    },
    # UKESM1-1 SSP2-4.5 from JASMIN. All three temperature variables sit in one
    # file, as UKESM G6-1.5K-SAI does, so the three entries point at the same
    # path and VAR_IN_FILE separates them. Coordinates are latitude/longitude
    # rather than lat/lon; see COORD_RENAME.
    ('UKESM', 'SSP245'): {
        'pr':     f'{_UKESM11}/{{ensemble_id}}/pr_day_UKESM1-1-LL_SSP245_{{ensemble_id}}_201501-210012.nc',
        'tas':    f'{_UKESM11}/{{ensemble_id}}/tas_mean_min_max_day_UKESM1-1-LL_SSP245_{{ensemble_id}}_201501-210012.nc',
        'tasmax': f'{_UKESM11}/{{ensemble_id}}/tas_mean_min_max_day_UKESM1-1-LL_SSP245_{{ensemble_id}}_201501-210012.nc',
        'tasmin': f'{_UKESM11}/{{ensemble_id}}/tas_mean_min_max_day_UKESM1-1-LL_SSP245_{{ensemble_id}}_201501-210012.nc',
    },

    # The UKESM1-0-LL reference already on GLADE. Unchanged paths, re-keyed
    # from 'UKESM' so the model version is in the key rather than in a note.
    ('UKESM10', 'SSP245'): {
        'pr':     f'{_UKESM}/SSP2-4.5/daily/pr/{{ensemble_id}}/',
        'tas':    f'{_UKESM}/SSP2-4.5/daily/tas/{{ensemble_id}}/',
        'tasmax': f'{_UKESM}/SSP2-4.5/daily/tasmax/{{ensemble_id}}/',
        'tasmin': f'{_UKESM}/SSP2-4.5/daily/tasmin/{{ensemble_id}}/',
    },

    # ---- E3SMv3 ------------------------------------------------------------
    ('E3SM', 'SAI'): None,
    ('E3SM', 'HiLLA'): None,
    ('E3SM', 'SSP245'): None,
}

# MODEL VERSION MISMATCH, UKESM. Confirmed 16 September 2026 from the CMOR
# global attributes of the GLADE files: UKESM SSP2-4.5 on GLADE is
# source_id 'UKESM1-0-LL' with variant labels r12i1p1f2, r2i1p1f2 and r3i1p1f2.
# The Reflective bucket's UKESM SSP2-4.5 is UKESM1-1 with r12i1p1f1 and
# r3i1p1f1: a different model version and a different forcing index. The
# G6-1.5K scenario filenames on GLADE say UKESM1-1-LL.
#
# So a UKESM build on GLADE would take its percentile threshold from one model
# version and apply it to another. That is not a path problem and no code
# change fixes it; it needs an answer from whoever staged the tree about which
# SSP2-4.5 the scenarios are meant to be scored against. Until then, treat any
# UKESM index built here as provisional.
VERSION_NOTE = {
    ('UKESM10', 'SSP245'): 'UKESM1-0-LL on GLADE (r12i1p1f2, r2i1p1f2, '
                           'r3i1p1f2); a different model version from the '
                           'G6-1.5K scenarios, kept for the version-sensitivity '
                           'comparison',
    ('UKESM', 'SSP245'): ('UKESM1.1 SSP2-4.5 from the JASMIN macloud '
                          'share: suites u-dp902, u-dp903 and u-dp904, '
                          'UM 13.1, stratospheric heterogeneous '
                          'chemistry on. The files carry no '
                          'variant_label or source_id; the labels r2, '
                          'r3 and r12 were confirmed from the UKESM '
                          'runs spreadsheet by Alistair Duffey on '
                          '25 September 2026.'),
}

# validation_note per (model, scenario), written by
# build_indices_derecho.py. Added 25 October 2026: the single sentence it
# used to stamp on every UKESM file cited different SSP2-4.5 runs in the
# two archives, untrue since 1 October 2026, and said nothing true about
# the G6 scenarios. repair_boot_attrs.py wrote the same sentences into
# the bootstrap trees.
VALIDATION_NOTE = {
    ('UKESM', 'SSP245'):
        ('The Reflective bucket holds copies of these SSP2-4.5 files '
         '(fafrifa/UKESM1-1/SSP245/, staged 1 October 2026 and verified by '
         'byte count), so this index can be checked against a Hub build '
         'member for member.'),
    ('UKESM', 'HiLLA'):
        ('The GLADE member directories 001, 002 and 003 hold suites '
         'u-ds287, u-ds275 and u-ds286 (variant labels r12i1p1f2, r2i1p1f2 '
         'and r3i1p1f2, read from the files on 1 October 2026), written here '
         'as r1, r2 and r3. The Hub reads CMORized files with the same '
         'variant labels. A Hub build of all 30 indices at 40-45N agreed with '
         'this tree to within float32 rounding (29-30 October 2026). The '
         'shared bucket copy of r3i1p1f2 stops before the end of 2084, so the '
         'GLADE copy was staged to the bucket (fafrifa/UKESM1-1/G6-1p5K-HiLLA/) '
         'on 30 October 2026.'),
    ('UKESM', 'SAI'):
        ('The GLADE member directories 001, 002 and 003 hold the '
         'stratospheric heterogeneous chemistry on suites u-dp691 '
         '(r12i1p1f2), u-dp583 (r2i1p1f2) and u-dp690 (r3i1p1f2), written '
         'here as r1, r2 and r3. They were identified on 26 October 2026 by '
         'matching daily tasmax against the suite-named files on the '
         'Reflective bucket (UKESM1-1/G6-1.5K-SAI/), and they pair with the '
         'het-chem-on SSP2-4.5 reference (u-dp902, u-dp903, u-dp904). The Hub '
         'reads the same suite files; a Hub build of all 30 indices at '
         '40-45N agreed with this tree to within float32 rounding '
         '(29 October 2026).'),
}

# No COORD_RENAME table. The JASMIN UKESM1-1 files store latitude/longitude
# rather than lat/lon, and one was added here on 19 September 2026 before the
# reader had been tested. preflight_ukesm11 then read SSP2-4.5 tas straight out
# of one of those files and got ('time', 'lat', 'lon') back, so derecho_io
# already renames this family of coordinates and the table did nothing. It is
# recorded as absent rather than silently dropped, because the next person to
# meet a latitude/longitude file will reach for exactly this and should know
# the renaming already happens one layer down.

# Internal variable names for the non-CMORized UKESM G6-1.5K-SAI files. The
# spreadsheet gives this mapping; the module already records that a UKESM member
# was excluded from the bucket for being uploaded non-CMORized with
# precipitation_flux and air_temperature_* names, so the same convention here is
# consistent, but the assignment of _0 and _1 is an assertion in a spreadsheet
# and nothing in a filename would reveal a swap. Verify against the file's own
# metadata before the first production run: read one file, print each variable's
# cell_methods and long_name, and confirm which carries maximum and which
# minimum. A swap of tasmax and tasmin produces TX90p and TN90p fields that look
# entirely plausible and are wrong.
VAR_IN_FILE = {
    ('UKESM', 'SAI'): {'tas': 'air_temperature_1',
                       'tasmax': 'air_temperature',
                       'tasmin': 'air_temperature_0'},
    # Same convention, and here it is not an assertion from a spreadsheet.
    # probe_ukesm11.py read cell_methods from all three files ('time: maximum',
    # 'time: minimum', 'time: mean') and checked it against the data: on the
    # first timestep, min <= mean <= max holds at every one of the 27,648 grid
    # points in each file. Metadata and values agree.
    ('UKESM', 'SSP245'): {'tas': 'air_temperature_1',
                          'tasmax': 'air_temperature',
                          'tasmin': 'air_temperature_0'},
}

# Unit conversions to mm/day, from the spreadsheet's notes.
PR_SCALE = {'CESM': 86400.0 * 1000.0,   # PRECT in m/s
            'MIROC': 86400.0,           # kg m-2 s-1
            'UKESM': 86400.0,
            'UKESM10': 86400.0,
            'E3SM': 86400.0}

# Land fraction files. MIROC and UKESM are percentages and need dividing by 100.
LAND_FRACTION = {
    'CESM': {'path': ('/glade/campaign/cesm/collections/ARISE-SAI-1.5/'
                      'b.e21.BW.f09_g17.SSP245-TSMLT-GAUSS-DEFAULT.001/atm/proc/'
                      'tseries/month_1/b.e21.BW.f09_g17.SSP245-TSMLT-GAUSS-'
                      'DEFAULT.001.cam.h0.LANDFRAC.203501-206912.nc'),
              'var': 'LANDFRAC', 'scale': 1.0, 'thresh': 0.5},
    'MIROC': {'path': f'{_MIROC}/sftlf_fx_*.nc',
              'var': 'sftlf', 'scale': 0.01, 'thresh': 0.5},
    # Note the file is UKESM1-0-LL piControl while the runs are UKESM1-1-LL.
    # Same N96 grid as far as I know, but confirm the lat/lon match before
    # trusting it, since a near-match would mask silently.
    'UKESM': {'path': f'{_UKESM}/sftlf_fx_UKESM1-0-LL_piControl_r1i1p1f2_gn.nc',
              'var': 'sftlf', 'scale': 0.01, 'thresh': 0.5},
}
# Both UKESM versions are N96, so they share one land fraction file.
LAND_FRACTION['UKESM10'] = LAND_FRACTION['UKESM']


def get(model, scenario, var):
    """Path template for one (model, scenario, variable), or raise with the
    reason it is missing."""
    key = (model, scenario)
    if key not in LAYOUT:
        raise KeyError(f'no Derecho layout entry for {key}')
    entry = LAYOUT[key]
    if entry is None:
        raise KeyError(
            f'{key} has no path in the team spreadsheet. Fill it from a real '
            f'listing on Derecho before running; do not derive it from a '
            f'neighbouring scenario.')
    if var not in entry:
        raise KeyError(f'{key} has no entry for variable {var!r}')
    return entry[var]


def members(model, scenario, through_year=None):
    """
    Member ids for a (model, scenario) pair.

    through_year drops members whose record ends before that year, which is how
    the two short MIROC HiLLA runs are excluded from the 2065-2084 window. Called
    without it, every member on disk is returned, which is what the baseline
    window wants.
    """
    key = (model, scenario)
    m = MEMBERS.get(key)
    if not m:
        raise KeyError(f'member ids for {key} are not recorded in the '
                       f'spreadsheet; list the directory to get them')
    if through_year is None:
        return m
    short = SHORT_MEMBERS.get(key, {})
    keep = tuple(x for x in m if short.get(x, 9999) >= through_year)
    dropped = [x for x in m if x not in keep]
    if dropped:
        print(f'  {key}: dropping {dropped} (record ends before {through_year})')
    return keep