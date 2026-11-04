# Selections in the `_tagsel.fits` file

**Status: proposal for review. Nothing here is implemented yet.** It replaces the single time window and single region
of selection script 2.x with any number of named time ranges and regions (selection script 3.0). Nobody depends on the
2.x layout, so there is no compatibility mode.

## In short

* Each time range or region owns one **bit** in a 16-bit flag column of the `EVENTS` table: `SEL_TIME` or `SEL_REGION`.
  An event inside several regions has several bits set; 0 means in none.
* Two table extensions, `TIMESEL` and `REGSEL`, define the bits: shape, frame, label, color, visibility.
* A region built from others (union, intersection, difference) is an **expression** over them and takes no bit.
* The existing `EDGE_OK` flag, positions, weights and processing keywords are unchanged.

## Flag columns

| Column | Type | Meaning |
|---|---|---|
| `SEL_TIME` | `I` with `TZERO = 32768` (unsigned 16-bit) | bit *n* set: the event is in time range with `BIT = n` |
| `SEL_REGION` | the same | bit *n* set: the event is in region with `BIT = n` |

Bit 0 is the least significant. astropy returns these as `uint16`; some other readers return signed values, in which
case bit 15 looks like a sign. At most 16 time ranges and 16 regions have a bit; derived rows do not count.
A column is written only if at least one selection of its kind exists.

## `TIMESEL` table (`EXTNAME = 'TIMESEL'`, `EXTVER = 1`)

| Column | Type | Meaning |
|---|---|---|
| `NAME` | string | short unique identifier (`T1`, `orbit2`): letters, digits, `_`; used in expressions |
| `NOTE` | string | the user's label, free text |
| `BIT`, `MASK` | `I`, `J` | bit number and its value `2**BIT`; `-1` and `0` for a derived row |
| `TSTART`, `TSTOP` | `D`, seconds | time range `TSTART <= TIME < TSTOP`, seconds since the exposure start (the `TIME` column); `NaN` for a derived row |
| `EXPR` | string | derived rows only: expression over earlier rows; empty otherwise |
| `EXPTIME` | `D`, seconds | good time (GTI) inside the range; derived, recomputable from the `GTI` extension |
| `COLOR` | string | `#rrggbb`, the color in plots |
| `VISIBLE` | `L` | shown in plots |

## `REGSEL` table (`EXTNAME = 'REGSEL'`, `EXTVER = 1`)

Same `NAME`, `NOTE`, `BIT`, `MASK`, `EXPR`, `COLOR`, `VISIBLE` columns, plus:

| Column | Type | Meaning |
|---|---|---|
| `FRAME` | string | reference frame the region was drawn in: `image`, `rest`, `rest_north`, `planet`, `planet_north` |
| `SHAPE` | string | `BOX`, `CIRCLE` or `POLYGON`; empty for a derived row |
| `PARAMS` | `nD`, pixel | shape numbers, in float64, padded with `NaN`; `n` is the longest row |
| `DESCR` | string | the one-line description the logs already use, for people |

`PARAMS` by shape:

| Shape | Numbers | Test |
|---|---|---|
| `BOX` | `x0, x1, y0, y1` | `x0 <= x < x1` and `y0 <= y < y1` |
| `CIRCLE` | `cx, cy, r` | `(x-cx)^2 + (y-cy)^2 <= r^2` |
| `POLYGON` | `x1, y1, x2, y2, ...` (closing edge implied) | inside by the even-odd rule; points exactly on an edge are not guaranteed either way |

A derived row has no `FRAME`; its frames are those of the rows in its expression. A new shape later (ellipse, annulus)
is a new `SHAPE` name with documented numbers; readers that do not know it can still use the bit column.

## Coordinates

* Regions use the **same 1-indexed hi-res pixel coordinates as the position columns** (`XCORR`, `XREST`, ...), so a
  region from the table can be applied by hand to those columns. The shift to 0-indexed arrays happens in one place in
  the package (the region is drawn at +1, plotted at −1, and one helper hands array-based code a 0-indexed copy).
* The flag is computed in float64 at each event's own position, before the position is rounded to the float32 column.
  Re-applying a region by hand to the float32 column can differ for the few events within about 1e-4 pixel of an edge.
  The flag is the authority.
* Planet-map frames (`map_*`, in degrees) are **not supported yet**: a region drawn there makes the script stop.

## Expressions

`EXPR` combines `NAME`s of rows **above it** in the same table, so there are no cycles:

| Operator | Meaning |
|---|---|
| `A \| B` | union |
| `A & B` | intersection |
| `A - B` | A and not B (same as `A & ~B`) |
| `~A` | not A |
| `( )` | grouping |

Precedence is `~`, then `&` and `-`, then `\|`; use parentheses to be clear. The text is parsed by the package, never
passed to `eval`. Time ranges and regions are separate tables and cannot refer to each other.
Combining across frames is exact, because each bit is already an event-level flag.

## What a script writes

* Position columns for the image frame (when corrections are on, or a region uses it) and for **every frame a region
  uses**, north-up or not, so each region can be re-applied by hand. `INCLUDE_NORTH_UP` still adds the rest.
* If any region's frame cannot be built (no SPT file, no Horizons answer), the script **stops before writing anything**.
* Header: `SELTMIN`, `SELTMAX`, `SELREG` and `SELFRAME` go away (the tables hold this); `TTETSEL` and `TTERSEL` stay,
  with the other `TTE*`, `RANDSEED`, `VIG*` and `SELEDGE*` keywords.
* `TIMESEL` and `REGSEL` follow the existing extensions, which are carried over untouched.

## Using it

The package will read the tables and hide the bit arithmetic:

```python
sel = read_selections('od8k57itq_tagsel.fits')            # time ranges and regions, with their numbers and colors
t = fits.getdata('od8k57itq_tagsel.fits', 1)
keep = in_time(t, 'orbit2') & in_region(t, 'oval') & t['EDGE_OK']   # by NAME; derived names work too
keep = (t['SEL_REGION'] & 0b101) != 0                      # or by hand: in the region with bit 0 or bit 2
```

Example tables for one exposure:

| NAME | NOTE | BIT | SHAPE / TSTART..TSTOP / EXPR |
|---|---|---|---|
| `orbit1` | first orbit | 0 | 0 to 1500 s |
| `orbit2` | second orbit | 1 | 3000 to 4500 s |
| `oval` | main oval | 0 | `POLYGON` in `planet_north` |
| `io` | Io footprint | 1 | `CIRCLE` in `rest` |
| `outside` | oval without footprint | −1 | `oval - io` |

"Brightness of the oval during orbit 2" is `in_time(t, 'orbit2') & in_region(t, 'oval')`. The file does not store such
pairings; the analysis picks them. (A `SELECTIONS` table that names them could be added later without changing any of
the above, for example to reload app state.)

## Open points

1. `NAME` is my addition so that expressions have something short and unique to refer to; `NOTE` is the free label.
2. 16 bits per column. A 32-bit word would remove the limit for twice the storage.
3. Map-frame regions (degrees) are left for a later version.
4. Spatial regions carry no exposure or area. Exposure varies across a region and depends on the time range and the
   corrections, so the package computes it on demand for a chosen pairing.
