import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from snapgene_parser import parse_snapgene_file
import os
import sys
import datetime


### Takes two dna sequences and compares their cut sites using the given enzyme library
# TODO - Currently only uses enzymes that cut in both sequences and ignores ones that don't cut in one.
   # This probably makes for a better screen as it confirms the vector is correct too so maybe no need to change

# TODO - make sequence a class

class Sequence:
    def __init__(self, filepath=None, sequence=None, name=None, circular=True):
        self._file_data = None
        self._circular = circular
        self._enzymes = {}
        if filepath:
            self._file_data = parse_snapgene_file(filepath)
            self._name = filepath[filepath.rfind('/') + 1:filepath.rfind('.')] if name is None else name
            self._circular = self._file_data["dna"]["topology"]
        else:
            if sequence is not None:
                self._file_data["seq"] = sequence
                self._name = name
            else:
                raise TypeError("Provide a string for at least one of filepath or sequence")

    def identify_enzymes(self):
        self._enzymes["cutters"], self._enzymes["non_cutters"] = find_enzyme_sites(self.sequence)

    @property
    def sequence(self):
        return self._file_data["seq"]

    @property
    def name(self):
        return self._name

    @property
    def circular(self):
        return self._circular

    def get_band_sizes(self, enzyme):
        if enzyme in self._enzymes["cutters"].keys():
            return get_band_sizes(self._enzymes["cutters"][enzyme], len(self.sequence))
        else:
            return [len(sequence)]



# Needed for the distributable to work
def get_aux_file_path(filename):
    if hasattr(sys, "_MEIPASS"):
        return os.path.join(sys._MEIPASS, "aux_files", filename)
    else:
        return os.path.join("aux_files", filename)


# TODO - This list is too basic, need cutting position as well for accurate sizing, as well as taking into account type II S enzymes (more than one cut per site)
elist_filename = get_aux_file_path("enzyme_list.csv")
ladders_filename = get_aux_file_path("ladders.csv")
enzyme_list = pd.read_csv(elist_filename)
ladders = pd.read_csv(ladders_filename, index_col=0)

sequence_validator = "ACGT"


def reverse_complement(seq):
    return seq[::-1].lower().replace('a', 'T').replace('t', 'A').replace('c', 'G').replace('g', 'C')


def validate_dna_sequence(seq):
    for char in seq:
        if char.upper() not in sequence_validator:
            print(char)
            return False
    return True


def sequence_comparison(sequence1: str, sequence2: str)->tuple:
    """
    Compares two sequences to identify potential enzymes that may cut differently
    :param sequence1: string of the first DNA sequence to compare
    :param sequence2: string of the second DNA sequence to compare
    returns a tuple of lists of enzymes that cut in just sequence 1, just sequence 2 or cut differently in sequence 1 and 2
    """
    sequence1_cutters, sequence1_non_cutters = find_enzyme_sites(sequence1)
    sequence2_cutters, sequence2_non_cutters = find_enzyme_sites(sequence2)
    cut1_not2 = []
    cut2_not1 = []
    cut_diff = []
    for k in sequence1_cutters.keys():
        if k in sequence2_non_cutters:
            cut1_not2.append({k: sorted(sequence1_cutters[k])})

        elif len(sequence1_cutters[k]) != len(sequence2_cutters[k]):
            cut_diff.append([(k, sorted(sequence1_cutters[k])), (k, sorted(sequence2_cutters[k]))])
    for k in sequence2_cutters.keys():
        if k in sequence1_non_cutters:
            cut2_not1.append({k: sorted(sequence2_cutters[k])})
    print(f"{len(cut1_not2)} enzymes found that cut in the first sequence, not the second")
    print(f"{len(cut2_not1)} enzymes found that cut in the second sequence, not the first")
    print(f"{len(cut_diff)} enzymes found that have a different number of cut sites in each sequence")
    return cut1_not2, cut2_not1, cut_diff
    ...


def get_band_sizes(site_list: list, plasmid_size: int=None)->list:
    """
    Gets the size of bands cut by an enzyme from a list of cut sites
    :param site_list: list of ints
    :param plasmid_size: if the sites refer to a plasmid then the plasmid size is needed to calculate correct band sizes
    returns a list of band sizes
    """
    band_sizes = []
    site_list.sort()
    if plasmid_size:
        b1 = site_list[0] + (plasmid_size - site_list[-1])
        band_sizes.append(b1)
    for i, site in enumerate(site_list, 1):
        if i == len(site_list):
            break
        band_sizes.append(site_list[i] - site)
    return sorted(band_sizes)


# For each band, if
def compare_band_sizes(band_list1: list, band_list2: list, cut_off:float=0.2)->list:
    """
    Compares two lists of integers to find the differences between the lists, use case is for comparing plasmids cut with the same enzymes to identify differences in the band length
    :param band_list1: list of integers
    :param band_list2: list of integers
    :param cut_off: the difference between the log10 values that should be returned
    returns a list of tuples of the index of band_list1 and band_list2 where differences are
    """
    min_band_diff = []
    log_b1 = np.log10(band_list1)
    log_b2 = np.vstack(np.log10(band_list2))
    band_diff = log_b1 - log_b2
    diffs_idx = np.where([np.all(i != 0) for i in band_diff])[0]
    diffs = np.abs(band_diff[diffs_idx])
    min_diffs_idx = [np.argmin(i) for i in diffs]
    assert len(diffs_idx) == len(min_diffs_idx)
    for i, _ in enumerate(diffs_idx):
        if np.abs(band_diff[diffs_idx[i]][min_diffs_idx[i]]) > cut_off:
            min_band_diff.append((diffs_idx[i], min_diffs_idx[i]))
    return min_band_diff


def find_enzyme_sites(sequence: str, circular: bool=True)->tuple:
    """
    Finds DNA sequence sites within a given sequence that are recognised by restriction enzymes
    :param sequence: string of dna sequence encompassing only A, G, C or Ts
    :param circular: set to True if the sequence is for a circular plasmid, set to False if the DNA sequence is linear
    returns a tuple of a dictionary of enzymes that cut and the cut site and a list of enzymes that don't cut
    """
    sequence = sequence.upper()
    rev_comp_sequence = reverse_complement(sequence)
    cutters = {}
    non_cutters = {}
    idx = enzyme_list.index
    for i in idx:
        e_name = enzyme_list.loc[i, "Enzyme"]
        e_seq = enzyme_list.loc[i, "Recognition Sequence"]
        match = [x for x in re.finditer(e_seq, sequence)]
        rev_match = [x for x in re.finditer(e_seq, rev_comp_sequence)]
        match2 = []
        rev_match2 = []
        if circular:
            match2 = [x for x in re.finditer(e_seq, sequence[-len(sequence) // 2:] + sequence[:-len(sequence) // 2])]
            rev_match2 = [x for x in re.finditer(e_seq, rev_comp_sequence[-len(rev_comp_sequence) // 2:] + rev_comp_sequence[:-len(rev_comp_sequence) // 2])]
        if len(match) > 0:
            cutters[e_name] = [x.start() for x in match]
        if len(rev_match) > 0:
            if e_name in cutters.keys():
                cutters[e_name] = list(set(cutters[e_name] + [(len(sequence) - x.end()) for x in rev_match]))
            else:
                cutters[e_name] = [len(sequence) - x.end() for x in rev_match]
        if len(match2) > 0:
            if e_name in cutters.keys():
                cutters[e_name] = list(set(cutters[e_name] + [(x.start() + (len(sequence) // 2)) % len(sequence) for x in match2]))
            else:
                cutters[e_name] = [(x.start() + (len(sequence) // 2)) % len(sequence) for x in match2]
        if len(rev_match2) > 0:
            if e_name in cutters.keys():
                cutters[e_name] = list(set(cutters[e_name] + [(len(sequence) - (x.end() + len(sequence) // 2)) % len(sequence) for x in rev_match2]))
            else:
                cutters[e_name] = [(len(sequence) - (x.end() + len(sequence) // 2)) % len(sequence) for x in rev_match2]
        if len(match) + len(rev_match) + len(match2) + len(rev_match2) == 0:
            non_cutters[e_name] = []
    return cutters, non_cutters


def terminal_gel_simulation(band_list):
    # band_list.sort(reverse=True)
    l = [x // 100 for x in band_list[::-1]]
    l  = [" " for _ in  range(100)]
    for x in band_list:
        l[x // 100] = "|"
    return "".join(l[::-1])
    # return '|'.join([' ' * x for x in l]) + "|"


def plot_gel(ladder, band_lists, enzyme):
    """
    Plots an image of the gel
    """
    fig, ax = plt.subplots()
    plot_gel_simulation(ladder, 0, ax)
    for i, band_list in enumerate(band_lists, 1):
        plot_gel_simulation(band_list, i, ax)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_xticks([0, *[i for i, _ in enumerate(band_lists, 1)]])
    ax.set_yscale("log")
    ax.set_yticks(ladder)
    ax.set_yticklabels(ladder)
    ax.set_ylim(100, 11000)
    ax.set_title(f"Digest with {enzyme}")
    return ax


def plot_gel_simulation(band_list, well, ax):
    ax.bar(
        [well for _ in band_list],
        [0 for _ in band_list],
        bottom=band_list,
        edgecolor="black",
        linewidth=5,
        )


def identify_best_cutter(sequence_1, sequence_2, enzyme_list):
    """
    Identifies the best options for a screening digest
    :param sequence_1: dna sequence for the first sequence
    :param sequence_2: dna sequence for the second sequence
    :param enzyme_list: list of enzymes given from the sequence_comparison function
    """
    best_cutters = []
    for n in enzyme_list:
        k1 = n[0][0]
        v1 = n[0][1]
        k2 = n[1][0]
        v2 = n[1][1]

        assert k1 == k2
        band_list1 = get_band_sizes(v1, plasmid_size=len(sequence_1))
        band_list2 = get_band_sizes(v2, plasmid_size=len(sequence_2))
        band_comparison = compare_band_sizes(band_list1, band_list2)
        if len(band_comparison) > 0:
            best_cutters.append((k1, band_list1, band_list2))
    return best_cutters


def get_dna_sequence_terminal(sequence_id: str=""):
    sequence = input(f"Enter {sequence_id} for screening digest identification: ")
    valid_sequence = validate_dna_sequence(sequence)
    while not valid_sequence:
        sequence = input("The sequence contains invalid characters.\n Please enter again just using A, C, G, T: ")
    return sequence


def get_sequence_from_snapgene(filepath):
    if filepath[0] in ["\"", "\'"]:
        filepath = filepath[1:]
    if filepath[-1] in ["\"", "\'"]:
        filepath = filepath[:-1]
    data = parse_snapgene_file(filepath)
    if data["dna"]["topology"] == "circular":
        return data["seq"]
    else:
        print("Currently only circular dna is supported")
        return None


def get_sequences_debug():
    return test_seqs["test_seq_1"], test_seqs["test_seq_2"]


def get_sequences():
    sequences = []
    add_sequence =  True
    if input("Do you want to use sequence from snapgene files (Y/N)? ") in "Yesyes":
        while add_sequence:
            temp_path = input("Enter full filepath to the sequence: ").replace("\'", "").replace("\"", "")
            sequences.append(Sequence(temp_path))
            # sequences.append(get_sequence_from_snapgene(temp_path))
            if input("Add more sequences (Y/N)? ") in "Nono":
                add_sequence = False
    else:
        while add_sequence:
            sequence_name = input("Enter name for the sequence: ")
            sequences.append(get_dna_sequence_terminal(sequence_name))
        sequence1 = Sequence(sequence=get_dna_sequence_terminal("Sequence 1"))
        sequence2 = Sequence(sequence=get_dna_sequence_terminal("Sequence 2"))
        sequences = [sequence1, sequence2]
    return sequences


def identify_cutters(sequences, reference_idx=0):
    """
    Identifys the enzymes that produce the best distinguishing cuts between sequences while cutting in both
    :params sequences: list of dna sequences
    :params reference_idx: gives the index of the reference sequence to which the others are compared
    """
    all_best_cutters = []
    # Identifys the best cutters for each sequence compared to a reference sequence
    for i, sequence in enumerate(sequences):
        if i == reference_idx:
            continue
        cut1, cut2, cut12 = sequence_comparison(sequences[reference_idx], sequence)
        best_cutters = [x[0] for x in identify_best_cutter(sequences[reference_idx], sequence, cut12)]

        # Identifying a common list between all sequences
        if (reference_idx == 0 and i == 1) or (i == 0):
            all_best_cutters = best_cutters
        else:
            temp_list = []
            for e in all_best_cutters:
                if e in best_cutters:
                    temp_list.append(e)
            all_best_cutters = temp_list
    return all_best_cutters


def get_bands_from_enzyme(enzyme, sequence, circular=True):
    """
    Given an enzyme and sequence, identify the cut sites
    """
    sequence = sequence.upper()
    rev_comp_sequence = reverse_complement(sequence)
    enzyme_loc = enzyme_list.loc[enzyme_list["Enzyme"] == enzyme]
    e_seq = enzyme_loc["Recognition Sequence"]
    cut_site_list = []

    match = [x for x in re.finditer(e_seq, sequence)]
    rev_match = [x for x in re.finditer(e_seq, rev_comp_sequence)]
    match2 = []
    rev_match2 = []

    if circular:
        match2 = [x for x in re.finditer(e_seq, sequence[-len(sequence) // 2:] + sequence[:-len(sequence) // 2])]
        rev_match2 = [x for x in re.finditer(e_seq, rev_comp_sequence[-len(rev_comp_sequence) // 2:] + rev_comp_sequence[:-len(rev_comp_sequence) // 2])]

    cut_site_list += [x.start() for x in match]
    cut_site_list += [len(sequence) - x.end() for x in rev_match]
    cut_site_list += [(x.start() + (len(sequence) // 2)) % len(sequence) for x in match2]
    cut_site_list += [(len(sequence) - (x.end() + len(sequence) // 2)) %  len(sequence) for x in rev_match2]
    return list(set(cut_site_list))


# '/Users/benjamindraper/Library/CloudStorage/Dropbox/01_UCL/Archive/Plasmids/[PL00309]pLKO.B7H3sgRNA1.mCherry.dna'
# '/Users/benjamindraper/Library/CloudStorage/Dropbox/01_UCL/Archive/Plasmids/[PL00279]pCCL[GGStuff].SFFV.RFP.dna'
# '/Users/benjamindraper/Library/CloudStorage/Dropbox/01_UCL/Archive/Plasmids/[PL00301]pCCL.TE9_BBz.dna'


def main():
    if debug:
        test_seqs = {}
        with open(get_aux_file_path("test.txt"), "r") as f:
            for i, line in enumerate(f.readlines(), 1):
                if len(line) > 1:
                    test_seqs[f"test_seq_{i}"] = line.strip()
        sequences = get_sequences_debug()
    else:
        sequences = get_sequences()

    for sequence in sequences:
        sequence.identify_enzymes()

    best_cutters = identify_cutters([seq.sequence for seq in sequences])

    print(f"{len(best_cutters)} Enzymes found that provide the best distinguishing bands: ")
    [print(f"[{i}] {cutter}") for i, cutter in enumerate(best_cutters)]
    ladder = list(ladders.loc["HyperLadderI"].values)
    show_gel = len(best_cutters) > 0
    while show_gel:
        show_gel_input = input("Enter the number of an enzyme to see a simulation of a gel, press q to quit: ")
        try:
            show_gel_int = int(show_gel_input)
            best_cutter = best_cutters[show_gel_int]

            bands = [seq.get_band_sizes(best_cutters[show_gel_int]) for seq in sequences]

            # Might be an error if get(best_cutter[0]) returns None, will need to have a way of showing that there is no cutting

            ax = plot_gel(ladder, bands, best_cutter)
            ax.set_xticklabels(["Ladder", *[seq.name for seq in sequences]], rotation=45)
            plt.show()


            # Can you save a matplotlib plot after plt.show()?
            if input("Save figure (Y/N)? ") in "Yesyes":
                save_path = input("Enter directory in which to save figure: ")
                if (len(save_path) > 0) and (save_path[-1] != "/"):
                    save_path = save_path + "/"
                date = datetime.datetime.now()
                iso_date = "-".join([format(date.year, "04"), format(date.month, "02"), format(date.day, "02")])
                ax = plot_gel(ladder, bands, best_cutter)
                ax.set_xticklabels(["Ladder", *[seq.name for seq in sequences]], rotation=45)
                plt.tight_layout()
                plt.savefig(f"{save_path}{iso_date}_{best_cutter}_digest.svg", bbox_inches="tight")
                plt.cla()
                print(f"Saved to {save_path} as {iso_date}_{best_cutter}_digest.svg")
            show_gel = input("Show another enzyme (Y/N)? ") in "Yesyes"
        except ValueError:
            if show_gel_input == "q":
                show_gel = False
            else:
                print("You made an error in selection, try again: ")


if __name__ == "__main__":
    debug = False
    main()





