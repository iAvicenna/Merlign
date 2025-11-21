#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 14:09:40 2023

@author: avicenna
"""
import numpy as np
import sys
import pandas as pd
import os 

file_path = os.path.dirname(os.path.realpath(__file__))

alph = ['A','G','C','T']
qalph = list('"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIZ') #! is reserved
ref0 = "GCTCAGGATGGGAAAAGCTATGCTTGCAAAAGGGGATCTGTTAACAGTTTCTTTAGTAGATTGAATTGGTTGCACAAATTAGAATACAAATATCCAGCGCTGAACGTGACTATGCCAAACAATGGCAAATTTGACAAATTGTACATTTGGGGGGTTCACCACCCGAGCACGGACAGTGACCAAACCAGCATATATGTTCGAGCATCAGGGAGAGTC"
MIDs = ["CATAGTAGTG", "ACGAGAGATAC", "CTACGTAGC", "TGTACTACTC", "CGACTACAG"]
adapters = ["AGATCGGAAGAG", "CTGTCTCTTATACAC"]

rc_map = {
  'A':'T',
  'C':'G',
  'T':'A',
  'G':'C',
  ' ':' '
  }

def reverse_complement(seq):

  return ''.join([rc_map[let] for let in seq[::-1]])

def generate_qual(length, rng):

  q = [qalph[rng.integers(0, len(qalph))] for _ in range(length)]

  return ''.join(q)


def mutation(seq, pos, rng, positions, subs_chance=0.8, ins_chance=0.01,
             del_chance=0, qual=None, _qual=None):

  dice = rng.random()
  nt0 = seq[pos]
  nt1 = list(set(alph).difference(nt0))[rng.integers(0,len(alph)-1)]

  if _qual is None:
    qt1 = qalph[rng.integers(0, len(qalph))]
  else:
    qt1 = _qual

  if dice>subs_chance + ins_chance + del_chance:
    return '', '', seq, qual, positions

  elif dice>subs_chance + ins_chance: #deletion

    if qual != None:
      qual = qual[:pos] + qual[pos+1:]

    positions = [x-1 if x>pos else x for x in positions]
    return ' ', ' ', seq[:pos] + seq[pos+1:], qual, positions

  elif dice>subs_chance: #insertion

    if qual != None:
      qual = qual[:pos] + qt1 + qual[pos:]

    positions = [x+1 if x>pos else x for x in positions]
    return f'{nt1}', f'{qt1}', seq[:pos] + nt1 + seq[pos:], qual, positions

  else: #substitution

    if qual != None:
      qual = qual[:pos] + qt1 + qual[pos+1:]

    return f'{nt1}', f'{qt1}', seq[:pos] + nt1 + seq[pos+1:], qual, positions


def generate_data(Nseqs, rng, nrefs=1, length=300):

  refs = [ref0]

  for i in range(nrefs-1):
    ref_new = ref0

    for j in range(200):
      _,_,ref_new,_,_ = mutation(ref_new, rng.integers(0, len(ref_new)), rng, [],
                                 ins_chance=0.1, del_chance=0.1)

    refs.append(ref_new)

  smods = ['' for _ in range(Nseqs)]
  qmods = ['' for _ in range(Nseqs)]
  ref_indices = []
  mid_indices = []
  adapter_indices = []
  is_rc = []

  reads=[]
  quals=[]
  ids=[]

  for i in range(Nseqs):

    ids.append(f'@Seq{i}')

    iref = rng.integers(0, len(refs))
    imid1 = rng.integers(0, len(MIDs))
    if rng.random()<0.2:
      imid2 = rng.integers(0, len(MIDs))
    else:
      imid2 = imid1
    iadap = rng.integers(0, len(adapters))

    ref = refs[iref]

    if rng.integers(0,2)==1:
      ref = reverse_complement(ref)
      is_rc.append(1)
      mid_indices.append((imid2, imid1))
    else:
      is_rc.append(0)
      mid_indices.append((imid1, imid2))

    mid1 = MIDs[imid1]
    mid2 = MIDs[imid2]
    adapter1 = adapters[iadap]
    adapter2 = adapter1
    q = generate_qual(len(ref), rng)

    ref_indices.append(iref)
    adapter_indices.append(iadap)

    offset = rng.integers(0,3)
    positions = rng.permutation(list(range(10+offset,len(ref)-10,3)))
    mut_positions = []

    smod, qmod, ref, q, positions = mutation(ref, positions[0], rng, positions,
                                             subs_chance=0.3, qual=q)

    if smod != '': mut_positions.append(positions[0])
    smods[i] += smod
    qmods[i] += qmod

    for j in range(1,min(51,len(positions))):
      smod, qmod, ref, q, positions = mutation(ref, positions[j], rng, positions,
                                               subs_chance=0.01, qual=q)

      if smod != '': mut_positions.append(positions[j])
      smods[i] += smod
      qmods[i] += qmod

    #since positions to be mutated are not in order above,
    #smods and qmods should be rearranged to have the correct order

    ordering = np.argsort(mut_positions)
    smods[i] = ''.join([smods[i][j] for j in ordering])
    qmods[i] = ''.join([qmods[i][j] for j in ordering])


    for j in range(2):
      qmid1 = generate_qual(len(mid1), rng)
      qmid2 = generate_qual(len(mid2), rng)
      qadap1 = generate_qual(len(adapter1), rng)
      qadap2 = generate_qual(len(adapter2), rng)

      _, _ , mid1, qmid1, _ = mutation(mid1, rng.integers(1, len(mid1)-1), rng, [],
                                       subs_chance=0.1, qual=qmid1)

      _, _ , mid2, qmid2, _ = mutation(mid2, rng.integers(1, len(mid2)-1), rng, [],
                                       subs_chance=0.1, qual=qmid2)

      _, _ , adapter1, qadap1, _ = mutation(adapter1, rng.integers(1, len(adapter1)-1),
                                          rng, [], subs_chance=0.1, qual=qadap1)

      _, _ , adapter2, qadap2, _ = mutation(adapter2, rng.integers(1, len(adapter2)-1),
                                          rng, [], subs_chance=0.1, qual=qadap2)


    read_forward = reverse_complement(adapter1) + mid1 + ref \
      + reverse_complement(mid2) + adapter2

    q_forward = qadap1[::-1] + qmid1 + q + qmid2 + qadap2

    reads.append(read_forward)
    quals.append(q_forward)

  read_records = [[id, read, qual] for id, read, qual in zip(ids, reads, quals)]

  return read_records, refs, smods, qmods, \
    [ref_indices, mid_indices, adapter_indices, is_rc]


def add_buffer(read_records, rng, length):

  new_read_records = []

  for id, read, qual in read_records:

    buffer = ''.join([alph[rng.integers(0, len(alph))] for _ in
                      range(0,length - len(read))])

    qbuf = generate_qual(len(buffer), rng)

    read = read + buffer
    qual = qual + qbuf

    new_read_records.append([id, read, qual])

  return new_read_records


def generate_pairs(read_records, is_rc, mid_indices, smods, qmods, rng):

  paired_records = [[],[],[]]

  offset = rng.integers(0,3)


  for indr, (read_id, read, qual) in enumerate(read_records):
    length = len(read)
    positions = rng.permutation(list(range(4+offset,length-10, 3)))

    new_qual = qual
    new_read = read
    for i in range(10):
      _, _, new_read, new_qual, [] = mutation(new_read, positions[i], rng, [],
                                              subs_chance=0.02, qual=new_qual,
                                               _qual='!', ins_chance=0)

    new_read = new_read[:length]
    new_qual = new_qual[:length]

    new_read = reverse_complement(new_read)
    new_qual = new_qual[::-1]

    paired_records[2].append([read_id, read, qual]) # only the read will be correct
                                                    # and I will only check the read
                                                    # in the test.

    if rng.integers(0,2)==0:
      i0=0;
      i1=1;
    else:
      i0=1;
      i1=0;

      is_rc[indr] += 1
      is_rc[indr] = is_rc[indr]%2

      smods[indr] = reverse_complement(smods[indr])
      qmods[indr] = qmods[indr][::-1]


    paired_records[i0].append([read_id, read, qual])
    paired_records[i1].append([read_id, new_read, new_qual])

  return paired_records


def write_records(records, path):

  with open(path, "w") as fp:
    for read_id, read, qual in records:
      fp.write(read_id+'\n')
      fp.write(read+'\n')

      if qual is not None:
        fp.write('+\n')
        fp.write(qual+'\n')

def create_output(indices, smods, qmods, read_ids):

  #[ref_indices, mid_indices, adapter_indices]

  mer_output = pd.DataFrame(-1, columns=["M_start", "M_end", "M_soft_start",	
                                         "M_soft_end", "M_match", "M_mismatch",
                                         "M_deletions", "M_insertions",
                                         "A1_index", "A1_start", "A1_end",
                                         "A1_match", "A1_mismatch", "A1_deletion",
                                         "A1_insertion", "A2_index", "A2_start", 
                                         "A2_end", "A2_match", "A2_mismatch", 
                                         "A2_deletion", "A2_insertion", "merged_len",
                                         "trimmed_len"], index=read_ids)

  ign_output = pd.DataFrame(-1, columns=["M1_index", "M1_start", "M1_end",
                                         "M1_match", "M1_mismatch", "M1_deletion", "M1_insertion",
                                         "M2_index", "M2_start", "M2_end", "M2_match",
                                         "M2_mismatch", "M2_deletion" , "M2_insertion",
                                         "align_index", "is_rc", "align_cigar",
                                         "align_mods", "align_mod_quals",
                                         "align_av_err"], index=read_ids);


  for indr,read_id in enumerate(read_ids):
    mer_output.loc[read_id, "A1_index"] = indices[2][indr]
    mer_output.loc[read_id, "A2_index"] = indices[2][indr]
    ign_output.loc[read_id, "M1_index"] = indices[1][indr][0]
    ign_output.loc[read_id, "M2_index"] = indices[1][indr][1]
    ign_output.loc[read_id, "align_index"] = indices[0][indr]
    ign_output.loc[read_id, "is_rc"] = indices[3][indr]
    ign_output.loc[read_id, "align_mods"] = smods[indr]
    ign_output.loc[read_id, "align_mod_quals"] = qmods[indr]


  return ign_output, mer_output


if __name__ == "__main__":

  if len(sys.argv) < 3:
    raise ValueError("Provide atleast two arguments, nseqs and nrefs (seed optional)")

  nseqs = int(sys.argv[1])
  nrefs = int(sys.argv[2])

  if len(sys.argv)==4:
    seed = int(sys.argv[3])
  else:
    seed = np.random.SeedSequence().spawn(1)[0]

  print(f"Generating {nseqs} sequences with {nrefs} references (seed={seed}).")

  rng = np.random.default_rng(seed)

  read_records, refs, smods, qmods, indices = generate_data(nseqs, rng, nrefs,
                                                            length=320)
  ref_records = [[f'>Ref{indr}',ref, None] for indr,ref in enumerate(refs) ]
  adapter_records = [[f'>Adapter{indr}',ref, None] for indr,ref in enumerate(adapters) ]
  is_rc = indices[-1]
  mid_indices = indices[1]
  paired_records = generate_pairs(read_records, is_rc, mid_indices, smods, qmods, rng)
  indices[-1] = is_rc
  indices[1] = mid_indices

  for i in range(2):
    paired_records[i] = add_buffer(paired_records[i], rng, 320)


  write_records(paired_records[0], f"{file_path}/../test_files/simulated_reads_R1.fastq")
  write_records(paired_records[1], f"{file_path}/../test_files/simulated_reads_R2.fastq")
  write_records(paired_records[2], f"{file_path}/../test_files/merged_reads.fastq")
  write_records(ref_records, f"{file_path}/../test_files/refs_for_simulated_reads.fasta")
  write_records(adapter_records, f"{file_path}/../test_files/adapters_for_simulated_reads.fasta")



  ign_output, mer_output = create_output(indices, smods, qmods,
                                         [x[0] for x in read_records])
  ign_output.to_csv(f"{file_path}/../test_files/simulated_ign_output.csv", index=True, header=True)
  mer_output.to_csv(f"{file_path}/../test_files/simulated_mer_output.csv", index=True, header=True)
