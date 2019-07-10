import itertools
import numpy as np
import os

from sequence_annotation_overlap import Sequences

path_to_script = os.path.dirname(os.path.abspath(__file__))
path_to_script = "/Users/elkeschaper/Documents/SIB/Tandem_repeats/stockholm_university/swissrepeat/src"

SWISSPROT_FILE = os.path.join(path_to_script, "../data/swissprot_annotations.tsv")
TR_ANNOTATIONS = os.path.join(path_to_script, "../results/tr_annotations/tr_annotations_post_filter.csv")
DISORDER_ANNOTATIONS = os.path.join(path_to_script, "../results/disorder_annotations/mobidb_regions_iupl.csv")
SEQUENCES = os.path.join(path_to_script, "../results/overview/annotation_data_trs_and_disorder.pickle")

if not os.path.exists(SEQUENCES):
    s = Sequences.create_sequences_from_csv(swissprot_filename=SWISSPROT_FILE)
    s.add_tr_annotations_from_csv(tr_filename=TR_ANNOTATIONS, tag="TR", id_name="ID")
    s.add_disorder_annotations_from_csv(disorder_filename=DISORDER_ANNOTATIONS, tag="IDR", delimiter=",")
    s.filter_annotations_with_minimum_overlap(tag_new="disordered_TR", tag_main="TR", tag_overlapping="IDR", overlap_ratio=0.9)
    s.filter_annotations_with_minimum_overlap(tag_new="disordered_TR_homo", tag_main="TR_homo", tag_overlapping="IDR", overlap_ratio=0.9)
    s.filter_annotations_with_minimum_overlap(tag_new="disordered_TR_micro", tag_main="TR_micro", tag_overlapping="IDR", overlap_ratio=0.9)
    s.filter_annotations_with_minimum_overlap(tag_new="disordered_TR_mini", tag_main="TR_mini", tag_overlapping="IDR", overlap_ratio=0.9)
    s.filter_annotations_with_minimum_overlap(tag_new="disordered_TR_long", tag_main="TR_long", tag_overlapping="IDR", overlap_ratio=0.9)
    s.write(file=SEQUENCES)
else:
    s = Sequences.create(file=SEQUENCES)


TAGS = ["TR", "TR_homo", "TR_micro", "TR_mini", "TR_long", "disordered_TR", "disordered_TR_homo", "disordered_TR_micro",
        "disordered_TR_mini", "disordered_TR_long", "IDR", "IDR_short", "IDR_medium", "IDR_long"]
GROUPS = ("superkingdom", ["Eukaryota", "Bacteria", "Archaea", "Viruses"])
STATISTICS = ["average_number_per_sequence", "number_of_sequences_with_tag", "relative_covered_chars"]

count_data = s.count_annotations(tags=TAGS, groupby_meta=GROUPS[0])
length_data = s.length_annotations(tags=TAGS, groupby_meta=GROUPS[0])
total_sequence_length_data = s.get_sequence_length(groupby_meta="superkingdom")

# Merge statistics :
merged_data = {tag: {i: {} for i in STATISTICS} for tag in TAGS}
for tag in TAGS:
    for group in GROUPS[1]:
        group_count_data = count_data[tag][group]
        merged_data[tag]["average_number_per_sequence"][group] = "{:.3f}".format(np.mean(group_count_data))
        merged_data[tag]["number_of_sequences_with_tag"][group]  = "{:.1f}".format(100 * np.count_nonzero(group_count_data)/len(group_count_data))
        #merged_data[tag]["relative_covered_chars"][group]  = "{:.2f}".format(100 * np.sum(length_data[tag]["coverage"][group])/total_sequence_length_data[group])
        merged_data[tag]["relative_covered_chars"][group]  = "{:.2f}".format(100 * np.mean([i for i in length_data[tag]["relative_coverage"][group] if i != 0]))


# Create table data:
tmp = [ "{}_{}".format(s, g) for s, g in itertools.product(STATISTICS, GROUPS[1])]
print(",".join(tmp))
for tag in TAGS:
    tmp = [merged_data[tag][s][g] for s, g in itertools.product(STATISTICS, GROUPS[1])]
    print(",".join(tmp))


# Interesting sequences for checking
special = []
for seq in s.sequences.values():
    if all([i in seq.annotations for i in ["TR_long", "IDR", "TR_mini", "TR_micro", "IDR_long"]]):
        special.append(seq.ID)
# special ['Q758C9', 'Q55FJ0', 'Q8S9G8', 'Q0A4N0', 'Q2HBI0', 'P52168', 'O75445', 'Q55FW7', 'Q3SYK4', 'Q55BP5', 'Q9V780', 'P12105', 'P9WIS7', 'B2HII5', 'Q8C1D8', 'Q3TLR7', 'O94966', 'Q62504', 'Q02447', 'P97479', 'O43451', 'Q02FS8', 'O62541', 'Q7Y0V7', 'Q86A14', 'Q3UTJ2', 'Q5RAA7', 'Q14517', 'Q3TY92', 'P9WNL8', 'Q9WTY8', 'Q6GGX3', 'Q5TM66', 'Q14678', 'Q6VMQ6', 'P05960', 'Q86L99', 'Q9NR30', 'Q5RAX9', 'P46097', 'O13493', 'Q9UBS5', 'Q3V0L5', 'Q5ZL34', 'Q21X09', 'B8AF63', 'Q5UQK6', 'P38536', 'Q9LFE4', 'P41180', 'Q07050', 'Q8K1Q4', 'A5VJE8', 'Q8NA61', 'Q54HP1', 'Q9NUG4', 'P0CM97', 'Q6Z1C0', 'Q9ZVD0', 'Q9VXE6', 'Q9AWK2', 'Q8WVC0', 'A1TLH9', 'Q99383', 'Q3L8U1', 'B4HEJ6', 'Q9SR13', 'P46889', 'P28166', 'B9W923', 'P14138', 'Q91VF5', 'Q96A84', 'P45386', 'Q0P5I0', 'O08852', 'Q8TD19', 'Q5AP95', 'Q6P4W5', 'Q08281', 'Q9LF24', 'Q7A5M1', 'Q9UKN7', 'D3ZKB9', 'Q79666', 'Q96MT7', 'Q19A40', 'Q9V3Q6', 'Q8IVF1', 'Q2FYJ6', 'Q8IUM7', 'Q292U2', 'Q9Z0P7', 'Q8CES0', 'O43281', 'A1A5S1', 'Q7NY13', 'Q4PSE8', 'P13383', 'Q86B20', 'Q54UL8', 'Q7T163', 'Q9P7E8', 'Q338N2', 'A9LNK9', 'P09493', 'Q96GP6', 'Q86HX1', 'Q6ZQF0', 'Q10951', 'Q5R439', 'Q55GT2', 'Q0ZM14', 'P04104', 'Q9JIT3', 'O08574', 'Q108T9', 'P32334', 'Q9LW88', 'Q99PJ1', 'Q2G0L4', 'Q53EP0', 'B0JU67', 'Q05470', 'D3KYU2', 'P22382', 'A8X3A7', 'Q3KSU8', 'E9Q9K5', 'P9WG29', 'Q6ZW49', 'Q5SP67', 'P0DKV0', 'P49824', 'Q96T23', 'Q1RMR2', 'Q5R684', 'Q9MAT0', 'A1KV51', 'Q9ZFC5', 'Q66JY2', 'Q04047', 'Q1L5Z9', 'Q54EQ9', 'P12524', 'Q4ZJI4', 'A0QRG5', 'Q501J7', 'F8RP11', 'Q9HSL6', 'Q9Z330', 'Q2J5Y1', 'Q6NXJ0', 'Q80U72', 'A1IGU3', 'P41391', 'O14335', 'Q6C8H2', 'Q8BLY3', 'Q9Y462', 'Q8SSW7', 'P34352', 'B2VR76', 'O75127', 'Q8BHE0', 'Q8N109', 'Q3V0C1', 'Q97T80', 'Q7M3S9', 'Q9Y6K5', 'Q2HAD8', 'Q921C3', 'Q9Z2L4', 'Q8N1F8', 'Q3UDK1', 'Q9NB71', 'Q5ALX5', 'P0CL52', 'Q9HCH0', 'P0C6K7', 'A4RN08', 'P0C7A2', 'A2A5K6', 'C4Y5P7']

for i in special[:10]:
    print(i)
    print(s.sequences[i].annotations)

