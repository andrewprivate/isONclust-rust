
extern crate parasailors;
extern crate needletail;
extern crate serde_json;
extern crate itertools;
use std::path::Path;
use std::fs;
use std::collections::HashMap;
use std::collections::VecDeque;
use parasailors::*;

// Will use traits later
pub struct InputSequence<'a> {
    read_cl_id: usize,
    batch_index: usize,
    sequence: &'a [u8],
    quality: &'a [u8],
    acc: &'a [u8], // accuracy?
    score: f64, // Score?
    error_rate: f64
}


pub fn cluster(sorted_data: &Vec<&InputSequence>) {
    // Arguments. Need to move out.
    let args_k = 1;
    let args_w = 1;
    let args_min_shared = 1;
    let args_min_shared = 1;
    let args_min_fraction = 1;
    let args_min_prob_no_hits = 1.0;
    let args_mapped_threshold = 1.0;

    // Will cache this later 
    let data = fs::read_to_string(Path::new("./p_minimizers_shared.json")).expect("Failed to read p_minimizers_shared file");
    let parsed: Vec<Vec<serde_json::Value>> = serde_json::from_str(&data).expect("Failed to parse p_minimizers_shared file");

   
    let mut p_emp_probs = HashMap::new();
    
    
    for values in &parsed {
        let k = values[0].as_i64().expect("Failed to get k");
        let w = values[1].as_i64().expect("Failed to get w");
        let p = values[2].as_f64().expect("Failed to get p");
        let e1 = values[3].as_f64().expect("Failed to get e1");
        let e2 = values[4].as_f64().expect("Failed to get e2");
        if k == args_k && (w - args_w).abs() <= 2 {
            // This is garbage, need to use proper data type. gonna use something like https://stackoverflow.com/questions/45786717/how-to-implement-hashmap-with-two-keys/45795699#45795699
            p_emp_probs.insert(e1.to_string() + "|" + &e2.to_string(), p);
            p_emp_probs.insert(e2.to_string() + "|" + &e1.to_string(), p);
        }
    }

    //  read_array = [ (i, 0, acc, seq, qual, float(acc.split("_")[-1])) for i, (acc, (seq, qual)) in enumerate(help_functions.readfq(open(sorted_reads_fastq_file, 'r')))]
    

    let mut clusters = HashMap::new();
    let mut representatives = HashMap::new();
    let mut minimizer_database = HashMap::new();

    for sequence in sorted_data {
        let mut n = Vec::new();
        n.push(sequence);
        clusters.insert(sequence.read_cl_id, n);
        let new = InputSequence { read_cl_id: sequence.read_cl_id, acc: sequence.acc, batch_index: sequence.batch_index, sequence: sequence.sequence, quality: sequence.quality, score: sequence.score, error_rate: 0.0 };
        representatives.insert(sequence.read_cl_id, new);
    }

    for sequence in sorted_data {

        // Step 1.1: Compress reads by removing consecutive duplicates
       // seq_hpol_comp = ''.join(ch for ch, _ in itertools.groupby(seq))
       
        //seq_hpol_comp.dedup();
        let mut all_read_hpol_lengths = Vec::new();
        let mut current_count = 0;
        let mut last_val: u8 = 0; // TODO, check if safe.
        let seq_hpol_comp: Vec<u8> = sequence.sequence.iter().filter(
            |val| {
                let is_duplicate = last_val == **val;
                last_val = **val;
                if is_duplicate {
                    current_count += 1;
                } else {
                    all_read_hpol_lengths.push(current_count);
                    current_count = 1;
                }
                return !is_duplicate;
            }
        ).cloned().collect();
        all_read_hpol_lengths.push(current_count);
       
        

        if seq_hpol_comp.len() < args_k as usize {
            print!( "skipping read of length: {} homopolymer compressed: {} sequence: {}", sequence.sequence.len(), seq_hpol_comp.len(), String::from_utf8(sequence.sequence.to_vec()).unwrap());
            continue;
        }

        // Step 1.2: Get minimizers
        let minimizers = get_kmer_minimizers(&seq_hpol_comp,args_k, args_w);

        // 2. Find the homopolymer compressed error rate (else statement is the only one active in single core mode)
       // if representatives.get(sequence.read_cl_id).unwrap().len() == 7 { // we have already computed homopolymenr compressed error rate in previous iteration (if isONclust is called with multiple cores):
            // not implemented as single-core only for now.
       //     panic!("Not implemented");
       // } else {
            let mut quality_comp = Vec::new();
            let mut start = 0;

            for h_len in all_read_hpol_lengths {
                let q_max = get_min_quality(&sequence.quality[start..(start + h_len)]);
                quality_comp.push(q_max);
                start += h_len;
            }
            assert_eq!(quality_comp.len(), seq_hpol_comp.len(), "Compressed quality length is not consistant with the compressed sequence length.");

            let quality_sum: f64 = quality_comp.iter().map(|v| v.1).sum();
            let h_pol_compr_error_rate = quality_sum/quality_comp.len() as f64;

        
            let new = InputSequence { read_cl_id: sequence.read_cl_id, acc: sequence.acc, batch_index: sequence.batch_index, sequence: sequence.sequence, quality: sequence.quality, score: sequence.score, error_rate: h_pol_compr_error_rate };
            representatives.insert(sequence.read_cl_id, new);
            
        //}
        // 3. Find all the representatives with shared minimizers (this is the time consuming function for noisy and large datasets)


        let (hit_clusters_ids, hit_clusters_hit_index, hit_clusters_hit_positions) = get_all_hits(&minimizers, &minimizer_database, sequence.read_cl_id);


        // # 4. Finds the best of the hits using mapping approach
        let (best_cluster_id_m, nr_shared_kmers_m, mapped_ratio) = get_best_cluster(sequence.read_cl_id, seq_hpol_comp.len(), &hit_clusters_ids, &hit_clusters_hit_positions, &minimizers, minimizers.len(), &hit_clusters_hit_index, &representatives, &p_emp_probs, args_min_shared, args_min_fraction, args_min_prob_no_hits, args_mapped_threshold);


        // 5. If step 4 is unsuccessfull we try to align the read to the representative(s) with the most shared minimizers.


        if best_cluster_id_m.is_none() && nr_shared_kmers_m >= args_min_shared {
            // aln_called += 1
            let (best_cluster_id_a, nr_shared_kmers_a, error_rate_sum, s1_alignment, s2_alignment, alignment_ratio) = get_best_cluster_block_align(sequence.read_cl_id, &representatives, &hit_clusters_ids, &hit_clusters_hit_positions);
            if best_cluster_id_a >= 0 {
              //  aln_passed_criteria += 1
            } else {
                best_cluster_id_a = -1
            }
        }
        
        // 6. Adds current read to representative, or makes it a new representative of a new cluster.
       



    }

}

fn get_best_cluster_block_align(read_cl_id: usize, representatives: &HashMap<usize, InputSequence>, hit_clusters_ids: &HashMap<usize, u64>, hit_clusters_hit_positions: &HashMap<usize, Vec<usize>>, args_k: i64) {
    let best_cluster_id = -1;
    let mut top_matches: Vec<(&usize, &Vec<usize>)> = hit_clusters_hit_positions.iter().collect();
    top_matches.sort_by_key(|x| (x.1.len(), x.1.iter().sum::<usize>(), representatives.get(x.0).unwrap().acc));
    top_matches.reverse();

    let rep_sequence = representatives.get(&read_cl_id).unwrap();
    let seq = rep_sequence.sequence;
    let r_qual = rep_sequence.quality;

    let top_hits = top_matches.get(0).unwrap().1.len();
    let alignment_ratio = 0.0;

    for tm in top_matches {
        let cl_id = tm.0;
        let nm_hits = tm.1.len();

        if nm_hits < top_hits {
            break;
        }

        let c_sequence = representatives.get(cl_id).unwrap();

        let c_seq = c_sequence.sequence;
        let c_qual = c_sequence.quality;
        
        let poisson_mean = sum_phred_char_to_p(r_qual);
        let poisson_mean2 = sum_phred_char_to_p(c_qual);

        let error_rate_sum = poisson_mean / seq.len() as f64 + poisson_mean2 / c_seq.len() as f64;

        let mut gap_opening_penalty = 0;
        if error_rate_sum <= 0.01 {
            gap_opening_penalty = 5;
        } else if error_rate_sum <= 0.04 {
            gap_opening_penalty = 4;
        } else if error_rate_sum <= 0.1 {
            gap_opening_penalty = 3;
        } else {
            gap_opening_penalty = 2;
        }

        let match_id_tailored = (1.0 - (error_rate_sum) * args_k as f64) as i64;
        
        let (s1, s2, (s1_alignment, s2_alignment, alignment_ratio)) = parasail_block_alignment(seq, c_seq, args_k, match_id_tailored, gap_opening_penalty);
      

    }
}

fn parasail_block_alignment(s1: &[u8], s2: &[u8], k: u64, match_id: i64, match_score: i64, mismatch_penalty: i64, opening_penalty: i64, gap_ext: i64) {

    // def parasail_block_alignment(s1, s2, k, match_id, match_score = 2, mismatch_penalty = -2, opening_penalty = 5, gap_ext = 1):

    let user_matrix = Matrix::create("ACGT", match_score, mismatch_penalty);
    let result = parasailors::

}

fn get_best_cluster(read_cl_id: usize, compressed_seq_len: usize, hit_clusters_ids: &HashMap<usize, u64>, hit_clusters_hit_positions: &HashMap<usize, Vec<usize>>, minimizers: &Vec<(&[u8], usize)>, nummber_of_minimizers: usize, hit_clusters_hit_index: &HashMap<usize, Vec<usize>>, representatives: &HashMap<usize, InputSequence>, p_emp_probs: &HashMap<String, f64>, args_min_shared: usize, args_min_fraction: usize, args_min_prob_no_hits: f64, args_mapped_threshold: f64) 
    -> (Option<usize>, usize, f64) {

    /*
          Tally up total covered (mapped region) and compare it with total unmapped region. What is counted as a consecutive mapped block depends on minimizer qualitity
        threshold as mapped-read-legnth/total-read-length is chosen to classify a read as belonging to the cluster
        
        Return: An integer >= 0 denoting the cluster ID that this read was assigned to. In not assigend to any previous cluster, return None. 
                [Also returing mapped ratio and nr shared minimizers to best read for logging purposes.]
    */
    let mut nr_shared_kmers = 0;
    let mut mapped_ratio = 0.0;
    let mut top_matches: Vec<(&usize, &Vec<usize>)> = hit_clusters_hit_positions.iter().collect();
    top_matches.sort_by_key(|x| (x.1.len(), x.1.iter().sum::<usize>(), representatives.get(x.0).unwrap().acc));
    top_matches.reverse();

    let top_hits = top_matches.get(0).unwrap().1.len();
    let nr_shared_kmers = top_hits;

    if top_hits < args_min_shared {
        // pass
    } else {
        for tm in top_matches {
            let cl_id = tm.0;
            let nm_hits = tm.1.len();
            if nm_hits < args_min_fraction * top_hits || nm_hits < args_min_shared {
                break;
            }

            let minimizer_hit_positions = hit_clusters_hit_positions.get(cl_id).unwrap();
            let minimizer_hit_indices = hit_clusters_hit_index.get(cl_id).unwrap();
            assert_eq!(minimizer_hit_indices.len(), minimizer_hit_positions.len(), "Unequal lengths!");

            let error_rate_c = representatives.get(cl_id).unwrap().error_rate;
            let error_rate_read = representatives.get(&read_cl_id).unwrap().error_rate;

            let p_error_in_kmers_emp = 1.0 - p_shared_minimizer_empirical(error_rate_read, error_rate_c, p_emp_probs);
            let minimizer_error_probabilities = vec![p_error_in_kmers_emp; nummber_of_minimizers];
            let total_mapped = 0;

            let prob_all_errors_since_last_hit = Vec::new();
            prob_all_errors_since_last_hit.push(minimizer_error_probabilities[..minimizer_hit_indices[0]].iter().fold(1.0, |a, b| a*b));
            for (hit_idx1, hit_idx2) in minimizer_hit_indices[..minimizer_hit_indices.len() - 1].iter().zip(minimizer_hit_indices[1..].iter()) {
                prob_all_errors_since_last_hit.push(minimizer_error_probabilities[hit_idx1+1..*hit_idx2].iter().fold(1.0, |a, b| a*b));
            }
            prob_all_errors_since_last_hit.push(minimizer_error_probabilities[minimizer_hit_indices[minimizer_hit_indices.len() - 1]+1..].iter().fold(1.0, |a, b| a*b));

            assert_eq!(prob_all_errors_since_last_hit.len(), minimizer_hit_positions.len() + 1, "Invalid lengths!");

            for i in 0..minimizer_hit_indices.len() {
                if prob_all_errors_since_last_hit[i] < args_min_prob_no_hits {
                    // pass
                } else {
                    if i == 0 {
                        total_mapped += minimizer_hit_positions[i];
                    } else {
                        total_mapped += minimizer_hit_positions[i] - minimizer_hit_positions[i - 1];
                    }
                }
            }
            if prob_all_errors_since_last_hit[prob_all_errors_since_last_hit.len() - 1] < args_min_prob_no_hits {
                // pass
            } else {
                total_mapped += compressed_seq_len - minimizer_hit_positions[minimizer_hit_positions.len() - 1];
            }

            mapped_ratio = total_mapped as f64 / compressed_seq_len as f64;

            if mapped_ratio > args_mapped_threshold {
                return (Some(*cl_id), nm_hits, mapped_ratio);
            }
            // [reduce(mul, minimizer_error_probabilities[: minimizer_hit_indices[0]], 1)] +  
            // [ reduce(mul, minimizer_error_probabilities[hit_idx1+1: hit_idx2], 1) for hit_idx1, hit_idx2 in zip(minimizer_hit_indices[:-1], minimizer_hit_indices[1:]) ] + 
            // [reduce(mul, minimizer_error_probabilities[minimizer_hit_indices[-1]+1 : ], 1)]

        }
    }
    
    return (None, nr_shared_kmers, mapped_ratio);
}

fn p_shared_minimizer_empirical(error_rate_read: f64, error_rate_center: f64, p_emp_probs: &HashMap<String, f64>) -> f64 {
    let mut e1 = (error_rate_read * 100.0).round() / 100.0;
    if e1 > 0.15 {
        e1 = 0.15;
    }
    if e1 < 0.01 {
        e1 = 0.01;
    }
    let mut e2 = (error_rate_center * 100.0).round() / 100.0;
    if e2 > 0.15 {
        e2 = 0.15;
    }
    if e2 < 0.01 {
        e2 = 0.01;
    }
    let p_kmer_shared = p_emp_probs.get(&(e1.to_string() + "|" + &e2.to_string())).unwrap();
    return *p_kmer_shared;
}

fn get_all_hits(minimizers: &Vec<(&[u8], usize)>, minimizer_database: &HashMap<&[u8], Vec<usize>>, read_cl_id: usize) 
    -> (HashMap<usize, u64>, HashMap<usize, Vec<usize>>,HashMap<usize, Vec<usize>>) {

    //     Get all representatives ID's that shares matches with the minimizers in the read.  

    let mut hit_clusters_ids = HashMap::new();
    let mut hit_clusters_hit_index = HashMap::new();
    let mut hit_clusters_hit_positions = HashMap::new();

    for (i, (m, pos)) in minimizers.iter().enumerate() {
        if minimizer_database.contains_key(m) {
            for cl_id in minimizer_database.get(m).unwrap() {
                hit_clusters_ids.insert(*cl_id, hit_clusters_ids.get(&cl_id).unwrap_or(&0) + 1);
                hit_clusters_hit_index.entry(*cl_id).or_insert(Vec::new()).push(i);
                hit_clusters_hit_positions.entry(*cl_id).or_insert(Vec::new()).push(*pos);
            }
        }
    }
    if hit_clusters_ids.contains_key(&read_cl_id) {
        hit_clusters_ids.remove(&read_cl_id);
        hit_clusters_hit_index.remove(&read_cl_id);
        hit_clusters_hit_positions.remove(&read_cl_id);
    }
    return (hit_clusters_ids, hit_clusters_hit_index, hit_clusters_hit_positions);
}

fn phred_char_to_p(c: u8) -> f64 {
    // phred_char_to_p = {chr(i) : min( 10**( - (ord(chr(i)) - 33)/10.0 ), 0.79433)  for i in range(128)} # PHRED encoded quality character to prob of error. Need this locally if multiprocessing
    return 10f64.powf(-(c as f64 - 33.0)/ 10.0).min(0.79433)
}
fn sum_phred_char_to_p(qual_str: &[u8]) -> f64 {
    return qual_str.iter().map(|x| phred_char_to_p(*x)).sum();
}
fn get_min_quality(qual: &[u8]) -> (u8, f64) {
  
    let mut min_val = qual[0];
    let mut min_sc = phred_char_to_p(min_val);
    for val in qual {
        let score = phred_char_to_p(*val);
        if score < min_sc {
            min_sc = score;
            min_val = *val;
        }
    }
    return (min_val, min_sc);
}

fn get_min_and_index<'a> (vec: &VecDeque<&'a [u8]>) -> (&'a [u8], usize) {
    
    let mut last_min = vec.get(0).unwrap();
    let mut last_index = 0;
    
    for (key, val) in vec.iter().enumerate() {
        if val < last_min {
            last_min = val;
            last_index = key;
        }
    }
    return (last_min, last_index);
}

fn get_kmer_minimizers(seq: &Vec<u8>, k_size: i64, w_size: i64) -> Vec<(&[u8], usize)> {
    let w = w_size - k_size;
    let mut window_kmers = VecDeque::new();

    // Make the kmers (sliding window)
    for i in 0..(w+1) {
        window_kmers.push_back(&seq[(i as usize)..((i+k_size) as usize)]);
    }

    let mut minimizers = Vec::new();
    let mut curr_min = get_min_and_index(&window_kmers);
    minimizers.push(curr_min);

    for i in (w as usize+1)..(seq.len() - k_size as usize + 1) {
        let new_kmer = &seq[i..(i + k_size as usize)];
        let discarded = window_kmers.pop_front().unwrap();
        window_kmers.push_back(new_kmer);

        if curr_min.0 == discarded {
            curr_min = get_min_and_index(&window_kmers);
            minimizers.push(curr_min);
        } else if new_kmer < curr_min.0 {
            curr_min = (new_kmer, i);
            minimizers.push(curr_min);
        }
    }
    return minimizers;
}


