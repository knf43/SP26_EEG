# noun acoustics
# Measures the duration, fundamental frequency (F0) and intensity of the noun
# in the two noun cue versions of each experimental item. The noun of
# conditions a, b and c is one recording (paroxytone cue) and the noun of
# conditions d, e and f is another (oxytone cue), so conditions a and d are
# measured. Times come from the trigger columns of the timing sheet.
#
# usage: python noun_acoustics.py <timing workbook> <sound folder> <output folder>

import os
import sys

import numpy as np
import openpyxl
import pandas as pd
import parselmouth
from scipy import stats

TIMING_SHEET = "5_experimental_psychopy_excel"
SENTENCE_SHEET = "2_experimental_sentences"
CUES = {"a": "paroxytone", "d": "oxytone"}
PITCH_FLOOR = 75
PITCH_CEILING = 400
TIME_STEP = 0.005
BIN = 0.050
LAST_BIN = 0.300
FINAL = 0.100


def read_sheet(workbook, sheet):
    rows = list(workbook[sheet].iter_rows(values_only=True))
    return pd.DataFrame(rows[1:], columns=rows[0])


def semitones(hertz):
    return 12 * np.log2(hertz / 100)


def window_mean(times, values, start, end):
    selected = values[(times >= start) & (times < end)]
    selected = selected[~np.isnan(selected)]
    return selected.mean() if selected.size >= 3 else np.nan


def measure(row, sound_folder):
    sound = parselmouth.Sound(os.path.join(sound_folder, os.path.basename(row["sounds"])))
    pitch = sound.to_pitch_ac(time_step=TIME_STEP, pitch_floor=PITCH_FLOOR,
                              pitch_ceiling=PITCH_CEILING)
    f0 = pitch.selected_array["frequency"]
    f0[f0 == 0] = np.nan
    f0_times = pitch.xs()
    intensity = sound.to_intensity(minimum_pitch=PITCH_FLOOR, time_step=TIME_STEP)
    db = intensity.values[0]
    db_times = intensity.xs()

    onset, f0_event, verb = (row["trigger_nounonset"], row["trigger_noun_F0"],
                             row["trigger_verbsonset"])
    in_noun = (f0_times >= f0_event) & (f0_times < verb) & ~np.isnan(f0)
    peak = np.argmax(np.where(in_noun, f0, -np.inf))

    result = {
        "item": row["item"], "noun": row["noun_s"], "cue": CUES[row["cond"]],
        "noun_duration": (verb - onset) * 1000,
        "onset_to_f0_event": (f0_event - onset) * 1000,
        "f0_event_to_verb": (verb - f0_event) * 1000,
        "f0_mean": semitones(window_mean(f0_times, f0, f0_event, verb)),
        "f0_peak": semitones(f0[peak]),
        "f0_peak_time": (f0_times[peak] - f0_event) * 1000,
        "f0_final": semitones(window_mean(f0_times, f0, verb - FINAL, verb)),
        "intensity_final": window_mean(db_times, db, verb - FINAL, verb),
    }
    for start in np.arange(0, LAST_BIN, BIN):
        label = f"{round(start * 1000)}_{round((start + BIN) * 1000)}"
        result[f"f0_{label}"] = semitones(
            window_mean(f0_times, f0, f0_event + start, f0_event + start + BIN))
        result[f"intensity_{label}"] = window_mean(
            db_times, db, f0_event + start, f0_event + start + BIN)
    return result


def main(timing_file, sound_folder, output_folder):
    workbook = openpyxl.load_workbook(timing_file, read_only=True, data_only=True)
    timing = read_sheet(workbook, TIMING_SHEET)
    sentences = read_sheet(workbook, SENTENCE_SHEET)[["sounds", "noun_s"]]
    timing = timing.merge(sentences, on="sounds")
    name = timing["sounds"].map(os.path.basename)
    timing["item"] = name.str.extract(r"x_(\d{2})")[0]
    timing["cond"] = name.str.extract(r"x_\d{2}([a-f])")[0]
    timing = timing[timing["cond"].isin(list(CUES))]

    items = pd.DataFrame([measure(row, sound_folder) for _, row in timing.iterrows()])
    items = items.round(2)
    items.to_csv(os.path.join(output_folder, "noun_acoustics_items.csv"), index=False)

    measures = [column for column in items.columns if column not in ("item", "noun", "cue")]
    wide = items.pivot(index="item", columns="cue", values=measures)
    rows = []
    for measure_name in measures:
        pair = wide[measure_name].dropna()
        difference = pair["oxytone"] - pair["paroxytone"]
        t = difference.mean() / (difference.std(ddof=1) / np.sqrt(len(difference)))
        rows.append({"measure": measure_name, "items": len(difference),
                     "paroxytone": pair["paroxytone"].mean(),
                     "oxytone": pair["oxytone"].mean(),
                     "difference": difference.mean(),
                     "sd_difference": difference.std(ddof=1), "t": t})
    summary = pd.DataFrame(rows)
    summary["p"] = 2 * stats.t.sf(summary["t"].abs(), summary["items"] - 1)
    summary = summary.round({"paroxytone": 2, "oxytone": 2, "difference": 2,
                             "sd_difference": 2, "t": 2, "p": 4})
    summary.to_csv(os.path.join(output_folder, "noun_acoustics_summary.csv"), index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main(*sys.argv[1:4])
