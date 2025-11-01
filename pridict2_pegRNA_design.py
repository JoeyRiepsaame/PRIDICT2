# -*- coding: utf-8 -*-
"""
Created on Tue May 11 16:04:15 2021

@author: m.wehling
"""

import pandas as pd
from Bio.Seq import Seq
import re
import time
import numpy as np
import os
import sys

# Get the absolute path of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Add the script's directory to the Python path
sys.path.append(script_dir)
from pridict_pegRNA_score import score_pegRNA

def primesequenceparsing(sequence):
    """Parse the input sequence and identify the desired edit."""
    original_base = re.search(r"\[([A-Z])", sequence)
    edited_base = re.search(r"\/([A-Z])", sequence)
    editposition_left = original_base.start()
    editposition_right = len(sequence) - original_base.start() - len(original_base.group())
    original_base = original_base.group(1)
    edited_base = edited_base.group(1)
    original_seq = sequence.replace('[', '').replace(']', '').replace('/' + edited_base, '')
    edited_seq = sequence.replace('[', '').replace(']', '').replace(original_base, edited_base, 1).replace('/' + edited_base,
                                                                                                           edited_base)
    mutation_type = original_base + '>' + edited_base
    correction_length = 1

    return original_base, edited_base, original_seq, edited_seq, editposition_left, editposition_right, mutation_type, correction_length, '', ''


def editorcharacteristics(editor):
    """Define characteristics of different prime editors."""
    if editor == 'PE2-NGG':
        PAM = '(?=GG)' # NGG PAM
        numberN = 2
        variant = 'SpCas9'
        protospacerlength = 20
        PAM_side = 'right'
        primescaffoldseq = 'GGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC'
        PAM_length = 3

    elif editor == 'PE2-NG':
        PAM = '(?=G)' # NG PAM
        numberN = 1
        variant = 'SpCas9-NG'
        protospacerlength = 20
        PAM_side = 'right'
        primescaffoldseq = 'GGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTTAACTTGCTATTTCTAGCTCTAAAAC'
        PAM_length = 2

    elif editor == 'PE2-SpRY':
        PAM = '(?=.)' # NRN PAM, but also works with NYN
        numberN = 1
        variant = 'SpCas9-SpRY'
        protospacerlength = 20
        PAM_side = 'right'
        primescaffoldseq = 'GGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC'
        PAM_length = 3
        
    elif editor == 'PE2-NRN':
        PAM = '(?=[AG].)' # NRN PAM (R = A or G)
        numberN = 1 
        variant = 'SpCas9-NRN'
        protospacerlength = 20
        PAM_side = 'right'
        primescaffoldseq = 'GGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC'
        PAM_length = 3

    return PAM, numberN, variant, protospacerlength, PAM_side, primescaffoldseq, PAM_length


def findpegRNAs(dfrow, PBSlengthrange=list(range(8, 18)), RToverhanglengthrange=list(range(1, 31)), windowsize_max=15, editor='PE2-NGG'):
    """Find pegRNAs and prediction scores for a set desired edit."""
    error_message = None
    try:
        if type(dfrow['editseq']) == pd.Series: # case of dfrow is a group
            sequence = dfrow['editseq'].values[0]
            name = dfrow['sequence_name'].values[0]
        else: # case of dfrow is a row in a dataframe
            sequence = dfrow['editseq'] # string
            name = dfrow['sequence_name']

        start_time = time.time()
        original_base, edited_base, original_seq, edited_seq, editposition_left, editposition_right, mutation_type, correction_length, basebefore_temp, baseafter_temp = primesequenceparsing(
            sequence)
        if (editposition_left < 99) or (editposition_right< 99):
            print('Less than 100bp flanking sequence! Check your input.')
            raise ValueError("Insufficient flanking sequence: Ensure at least 100 bp flanking the edit on both sides.")
            
        sequence = sequence.upper()
        PAM, numberN, variant, protospacerlength, PAM_side, primescaffoldseq, PAM_length = editorcharacteristics(editor)

        mutationtypelist = []
        correctiontypelist = []
        correctionlengthlist = []
        edited_sequence_list = []
        revcomp_edited_sequence_list = []
        original_sequence_list = []
        revcomp_original_sequence_list = []
        mutation_position_to_PAM = []
        editedallelelist = []
        originalallelelist = []
        variantList = []
        target_strandList = []
        protospacerpamsequence = []
        protospacer_oligo_FW = []
        protospacer_oligo_RV = []
        extension_oligo_FW = []
        extension_oligo_RV = []
        pcrGG_oligo1 = []
        pcrGG_oligo2 = []
        editpositionlist = []
        PBSlength_variants_dic = {}
        PBSrevcomp_dic = {}
        for length in PBSlength_variants:
            PBSlength_variants_dic[length] = []
            PBSrevcomp_dic[length] = []
        RTseqoverhang_variants_dic = {}
        for length in RTseqoverhang_variants:
            RTseqoverhang_variants_dic[length] = []

        PBSsequencelist = []
        PBSrevcomplist = []
        RTseqlist = []
        RTseqoverhangrevcomplist = []
        RTseqrevcomplist = []
        deepcas9seqlist = []
        pbslengthlist = []
        rtlengthoverhanglist = []
        rtlengthlist = []
        pegRNA_list = []
        scaffold_list = []
        nickingprotospacerlist = []
        nickingpositiontoeditlist = []
        nickingprotospacerlist = []
        nickingdeepcas9list = []
        nickingpe3blist = []
        nickingPAMdisruptlist = []
        nickingtargetstrandlist = []
        nicking_oligo_FW = []
        nicking_oligo_RV = []
        Doench2016_Score_list = []
        pegRNA_score_list = []
        target_sequence_list = []

        for target_strand, windowsize in {'FW': windowsize_max, 'RV': windowsize_max}.items():
            editposition = editposition_left
            if target_strand == 'RV':
                original_seq = str(Seq(original_seq).reverse_complement())
                edited_seq = str(Seq(edited_seq).reverse_complement())

            editingWindow = range(0 - windowsize + 3, 4)
            temp_dic = {}

            X = [m.start() for m in re.finditer(PAM, original_seq)]
            X = [x for x in X if 25 <= x < len(original_seq) - 4]
            
            editedstrand_PAMlist = [m.start() for m in re.finditer(PAM, edited_seq.upper())]
            editedstrand_PAMlist = [x for x in editedstrand_PAMlist if 25 <= x < len(edited_seq) - 4]

            if X:
                xindex = 0
                editedstrandPAMindex = 0
                for editedstrandPAM in editedstrand_PAMlist:
                    editedPAM_int = editedstrand_PAMlist[editedstrandPAMindex] - editposition - numberN
                    editedstrandPAMindex = editedstrandPAMindex + 1
                    editedPAM = editedPAM_int + editposition + numberN
                    editedstart = editedPAM + (len(PAM) - 7) - 3
                    editednickposition = editedstart - editposition
                    editednickprotospacer, editednickdeepcas9 = nickingguide(edited_seq, editedstrandPAM, protospacerlength)
                    nickingtargetstrandlist.append(target_strand)
                    nickingpositiontoeditlist.append(editednickposition)
                    nickingprotospacerlist.append(editednickprotospacer)
                    nickingdeepcas9list.append(editednickdeepcas9)
                    if editednickprotospacer[0] != 'G':
                        editednickprotospacer = 'g'+editednickprotospacer
                    nicking_oligo_FW.append('cacc' + editednickprotospacer)
                    nicking_oligo_RV.append('aaac' + str(Seq(editednickprotospacer).reverse_complement()))
                    pe3bwindowlist = list(range(-5, 17))
                    pe3bwindowlist.remove(-3)

                    if editednickposition in pe3bwindowlist:
                        nickingpe3blist.append('PE3b')
                        if editednickposition in [-5,-4]:
                            nickingPAMdisruptlist.append('Nicking_PAM_disrupt')
                        else:
                            nickingPAMdisruptlist.append('No_nicking_PAM_disrupt')
                    else:
                        nickingpe3blist.append('No_PE3b')
                        nickingPAMdisruptlist.append('No_nicking_PAM_disrupt')

                for xvalues in X:
                    X_int = X[xindex] - editposition - numberN
                    xindex = xindex + 1
                    XPAM = X_int + editposition + numberN
                    start = XPAM + (len(PAM) - 7) - 3
                    if X_int in editingWindow:
                        protospacer, deepcas9seq = spacer(original_seq, XPAM, protospacerlength)
                        for PBS_length in PBSlength_variants:
                            PBS, PBSrevcomp = pbs(edited_seq, start, PBS_length)
                            PBSlength_variants_dic[PBS_length].append(PBS)
                            PBSrevcomp_dic[PBS_length].append(PBSrevcomp)
                        for RT_overhang in RTseqoverhang_variants:
                            RTseq, RTseqrevcomp = rt(edited_seq, start, correction_length, RT_overhang)
                            RTseqoverhang_variants_dic[RT_overhang].append(RTseq)
                        
                        Doench2016_Score, pegRNA_score, target_sequence = score_pegRNA(protospacer, deepcas9seq, RTseq, PBS)
                        Doench2016_Score_list.append(Doench2016_Score)
                        pegRNA_score_list.append(pegRNA_score)
                        target_sequence_list.append(target_sequence)

                        mutationtypelist.append(mutation_type)
                        correctionlengthlist.append(correction_length)
                        edited_sequence_list.append(edited_seq)
                        original_sequence_list.append(original_seq)
                        mutation_position_to_PAM.append(X_int)
                        editedallelelist.append(edited_base)
                        originalallelelist.append(original_base)
                        variantList.append(variant)
                        target_strandList.append(target_strand)
                        protospacerpamsequence.append(protospacer)
                        editpositionlist.append(editposition)

                        if protospacer[0] != 'G':
                            protospacer = 'g' + protospacer
                        protospacer_oligo_FW.append('cacc' + protospacer)
                        protospacer_oligo_RV.append('aaac' + str(Seq(protospacer).reverse_complement()))
                        extension_oligo_FW.append('g' + RTseqrevcomp + PBSrevcomp)
                        extension_oligo_RV.append(str(Seq('g' + RTseqrevcomp + PBSrevcomp).reverse_complement()))
                        pcrGG_oligo1.append('CCTTCTAGAA' + protospacer)
                        pcrGG_oligo2.append('GTTTAAGAGC' + str(Seq(PBSrevcomp).reverse_complement()) + str(
                            Seq(RTseqrevcomp).reverse_complement()))

        df = pd.DataFrame()
        if not protospacerpamsequence:
            error_message = 'No PAM found in proximity of the edit'
            print('No PAM found in proximity of the edit')
        else:
            df['pegRNA sequence'] = protospacerpamsequence
            df['pegRNA score'] = pegRNA_score_list
            df['Doench2016 score'] = Doench2016_Score_list
            df['Target sequence'] = target_sequence_list
            df['Mutation type'] = mutationtypelist
            df['Correction length'] = correctionlengthlist
            df['Position of PAM to edit'] = mutation_position_to_PAM
            df['Cas-variant'] = variantList
            df['Target strand'] = target_strandList
            df['Protospacer_oligo_FW'] = protospacer_oligo_FW
            df['Protospacer_oligo_RV'] = protospacer_oligo_RV
            df['Extension_oligo_FW'] = extension_oligo_FW
            df['Extension_oligo_RV'] = extension_oligo_RV
            df['original_sequence'] = original_sequence_list
            df['edited_sequence'] = edited_sequence_list
            df['edit_position'] = editpositionlist
            
            for pbslength in PBSlength_variants:
                df['PBS_length_' + str(pbslength)] = PBSlength_variants_dic[pbslength]
            for rtlength in RTseqoverhang_variants:
                df['RT_length_' + str(rtlength)] = RTseqoverhang_variants_dic[rtlength]

            df2 = pd.DataFrame()
            df2['Nicking_guide_protospacer'] = nickingprotospacerlist
            df2['Nicking_guide_position'] = nickingpositiontoeditlist
            df2['PE3b'] = nickingpe3blist
            df2['Nicking_guide_PAM_disruption'] = nickingPAMdisruptlist
            df2['Nicking_guide_strand'] = nickingtargetstrandlist
            df2['Nicking_guide_oligo_FW'] = nicking_oligo_FW
            df2['Nicking_guide_oligo_RV'] = nicking_oligo_RV

            df['sequence_name'] = name
            df2['sequence_name'] = name
            df = df.sort_values(by='pegRNA score', ascending=False)
            df2 = df2.sort_values(by='Nicking_guide_position', ascending=True)

        return df, df2, error_message

    except Exception as e:
        error_message = f"An error occurred: {e}"
        return pd.DataFrame(), pd.DataFrame(), error_message


def spacer(sequence, x, protospacerlength):
    """Get protospacer and DeepSpCas9 input sequence."""
    protospacer = sequence[x - protospacerlength:x]
    deepcas9seq = sequence[x - protospacerlength - 4:x + 3 + 3]
    return protospacer, deepcas9seq


def nickingguide(sequence, x, protospacerlength):
    """Get nicking guide and DeepSpCas9 input sequence."""
    protospacer = sequence[x - protospacerlength:x]
    deepcas9seq = sequence[x - protospacerlength - 4:x + 3 + 3]
    return protospacer, deepcas9seq


def pbs(sequence, start, length):
    """Get PBS sequence."""
    PBS = sequence[start - length:start]
    PBSrevcomp = str(Seq(PBS).reverse_complement())
    return PBS, PBSrevcomp


def rt(sequence, start, correction_length, overhang):
    """Get RT sequence."""
    RTseq = sequence[start:start + correction_length + overhang]
    RTseqrevcomp = str(Seq(RTseq).reverse_complement())
    return RTseq, RTseqrevcomp
