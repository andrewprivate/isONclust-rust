
extern crate parasailors;
extern crate needletail;
extern crate serde_json;
extern crate itertools;
extern crate regex;
extern crate ordered_float;
use regex::Regex;
use std::path::Path;
use std::fs;
use std::collections::HashMap;
use std::collections::VecDeque;
use parasailors::*;
use std::cmp;

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

pub fn process_fastq(reads: Vec<(String, String, String)>, k: usize) {
    let mut probs = HashMap::new();
    for i in 33u8..128 {
        probs.insert(i, 10f64.powf(-1f64 * ((i as f64 - 33.0) / 10.0)));
        // println!("{}\t{}", i as char, 10f64.powf(-1f64 * ((i as f64 - 33.0) / 10.0)));
    }

    let mut sorted_data = Vec::new();
    for (i, (acc, sequence, qual)) in reads.iter().enumerate() {
        if sequence.len() < 2 * k {
            continue
        }
        // println!("{}", acc);
        let score = calculate_score(&qual, k, &probs);
        // println!("Score: {}", score);
        let mut new = InputSequence{ read_cl_id: i, batch_index: 0, sequence: sequence.as_bytes(), quality: qual.as_bytes(), acc: acc.as_bytes(), score, error_rate: 0f64};
        sorted_data.push(new);
    }
    sorted_data.sort_by(|s1, s2| s2.score.partial_cmp(&s1.score).unwrap());
    // for r in sorted_data.iter() {
        // println!("@{}_{}\n{}\n+\n{}", String::from_utf8_lossy(r.acc), r.score, String::from_utf8_lossy(r.sequence), String::from_utf8_lossy(r.quality));
        // println!("@{}_{}", String::from_utf8_lossy(r.acc), r.score);
    // }
    cluster(&sorted_data);
}

fn calculate_score(qual: &str, k: usize, probs: &HashMap<u8, f64>) -> f64 {
    // println!("Len: {}", qual.len());
    let exp_errors_in_kmers = estimate_erroneous_kmers(qual, k, probs);
    // println!("Exp errors: {}", exp_errors_in_kmers);
    let p_no_errors = 1f64-(exp_errors_in_kmers)/(qual.len() - k + 1) as f64;
    // println!("P no errors: {}", p_no_errors);
    let score = p_no_errors * (qual.len() - k + 1) as f64;

    return score;
}

fn estimate_erroneous_kmers(qual: &str, k: usize, probs: &HashMap<u8, f64>) -> f64 {
    let mut prob_vec = Vec::new();
    for c in qual.chars() {
        prob_vec.push(1. - *probs.get(&(c as u8)).unwrap());
    }
    // println!("{:?}", prob_vec);
    let mut window: VecDeque<f64> = prob_vec[0..k].iter().cloned().collect();
    // println!("{:?}", window);
    let mut prob_no_error = 1f64;
    for p in &window {
        prob_no_error *= *p;
    }
    // println!("{}", prob_no_error);
    let mut sum_of_e = prob_no_error;
    for p in prob_vec[k..].iter() {
        let p_to_leave = window.pop_front().unwrap();
        prob_no_error *= p/p_to_leave;
        sum_of_e += prob_no_error;
        window.push_back(*p);
    }
    return (qual.len() - k + 1) as f64 - sum_of_e
}

pub fn build_input_sequence<'a> ( read_cl_id: usize, batch_index: usize, sequence: &'a [u8], quality: &'a [u8], acc: &'a [u8],  score: f64, error_rate: f64) -> InputSequence<'a> {
    InputSequence {
        read_cl_id,
        batch_index,
        sequence,
        quality,
        acc, // accuracy?
        score, // Score?
        error_rate
    }
}

pub fn cluster<'a>(sorted_data: &'a Vec<InputSequence<'a>>) -> (HashMap<usize, std::vec::Vec<&InputSequence<'a>>>, HashMap<usize, InputSequence<'a>>, HashMap<Vec<u8>, std::vec::Vec<usize>>) {
    // Arguments. Need to move out.
    let args_k = 15;
    let args_w = 50;
    // let args_min_shared = 5;
    let args_min_shared = 5;
    let args_min_fraction = 0.8;
    let args_aligned_threshold = 0.4;
    let args_min_prob_no_hits = 0.1;
    let args_mapped_threshold = 0.7;
    let mut aln_calls = 0;

    // Will cache this later 
    let data = fs::read_to_string(Path::new("C:\\CLionProjects\\isONclust-rust\\src\\p_minimizers_shared.json")).expect("Failed to read p_minimizers_shared file");
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

    let mut cluster_to_new_cluster_id = HashMap::new();

    for sequence in sorted_data {

        // Step 1.1: Compress reads by removing consecutive duplicates
       // seq_hpol_comp = ''.join(ch for ch, _ in itertools.groupby(seq))
       
        //seq_hpol_comp.dedup();
        
        let (seq_hpol_comp, all_read_hpol_lengths) = get_seq_hpol_comp(sequence.sequence);
        if seq_hpol_comp.len() < args_k as usize {
            eprintln!( "skipping read of length: {} homopolymer compressed: {} sequence: {}", sequence.sequence.len(), seq_hpol_comp.len(), String::from_utf8(sequence.sequence.to_vec()).unwrap());
            continue;
        }

        // Step 1.2: Get minimizers
        let minimizers = get_kmer_minimizers(&seq_hpol_comp, args_k, args_w);
        // 2. Find the homopolymer compressed error rate (else statement is the only one active in single core mode)
       // if representatives.get(sequence.read_cl_id).unwrap().len() == 7 { // we have already computed homopolymenr compressed error rate in previous iteration (if isONclust is called with multiple cores):
            // not implemented as single-core only for now.
       //     panic!("Not implemented");
       // } else {
            let mut quality_comp = Vec::new();
            let mut start = 0;

            for h_len in all_read_hpol_lengths {
                if h_len > 0 {
                    let q_max = get_min_quality(&sequence.quality[start..(start + h_len)]);
                    quality_comp.push(q_max);
                    start += h_len;
                }
            }
            assert_eq!(quality_comp.len(), seq_hpol_comp.len(), "Compressed quality length is not consistant with the compressed sequence length.");

            let quality_sum: f64 = quality_comp.iter().map(|v| v.1).sum();
            let h_pol_compr_error_rate = quality_sum/quality_comp.len() as f64;
            // println!("{}", h_pol_compr_error_rate);
        
            let new = InputSequence { read_cl_id: sequence.read_cl_id, acc: sequence.acc, batch_index: sequence.batch_index, sequence: sequence.sequence, quality: sequence.quality, score: sequence.score, error_rate: h_pol_compr_error_rate };
            representatives.insert(sequence.read_cl_id, new);

        //}
        // 3. Find all the representatives with shared minimizers (this is the time consuming function for noisy and large datasets)


        let (hit_clusters_ids, hit_clusters_hit_index, hit_clusters_hit_positions) = get_all_hits(&minimizers, &minimizer_database, sequence.read_cl_id);


        // # 4. Finds the best of the hits using mapping approach
        let (best_cluster_id_m, nr_shared_kmers_m, _mapped_ratio) = get_best_cluster(sequence.read_cl_id, seq_hpol_comp.len(), &hit_clusters_ids, &hit_clusters_hit_positions, &minimizers, minimizers.len(), &hit_clusters_hit_index, &representatives, &p_emp_probs, args_min_shared, args_min_fraction, args_min_prob_no_hits, args_mapped_threshold);


        // 5. If step 4 is unsuccessfull we try to align the read to the representative(s) with the most shared minimizers.


        let mut best_cluster_id_a = -1;

        if best_cluster_id_m < 0 && nr_shared_kmers_m >= args_min_shared {
            aln_calls += 1;
            let (best_cluster_id_ta, _nr_shared_kmers_a, _error_rate_sum, _s1_alignment, _s2_alignment, _alignment_ratio) = get_best_cluster_block_align(sequence.read_cl_id, &representatives, &hit_clusters_ids, &hit_clusters_hit_positions, args_k, args_aligned_threshold);
            best_cluster_id_a = best_cluster_id_ta;
        }
        
        // 6. Adds current read to representative, or makes it a new representative of a new cluster.
       
        let best_cluster_id = cmp::max(best_cluster_id_m, best_cluster_id_a);
       
        if best_cluster_id >= 0 {
            cluster_to_new_cluster_id.insert(sequence.read_cl_id, best_cluster_id as usize);
        } else {
            // 7. If new representative: add the minimizers to the minimizer database. 

            for (k, pos) in minimizers {
                if minimizer_database.contains_key(k) {
                    let md_vec = minimizer_database.get_mut(k).unwrap();
                    if md_vec.contains(&sequence.read_cl_id) {

                    } else {
                        md_vec.push(sequence.read_cl_id);
                    }
                } else {
                    let mut v = Vec::new();
                    v.push(sequence.read_cl_id);
                    minimizer_database.insert(k.to_owned(), v);
                }
            }
        }
    }

    for (old_cl_id,new_cl_id) in cluster_to_new_cluster_id {
        // let new_cl_id = cluster_to_new_cluster_id.get(read_cl_id).unwrap();
        let all_reads = clusters.get(&old_cl_id).unwrap().to_owned();
        let cl = clusters.get_mut(&new_cl_id).unwrap();
        for read_acc in all_reads {
            cl.push(read_acc);
        }
        clusters.remove(&old_cl_id);
        representatives.remove(&old_cl_id);
    }

    let mut total = 0;
    let mut g1 = 0;
    for (cl_id, read_vec) in &clusters {
        total += 1;
        if read_vec.len() > 1 {
            g1 += 1;
        }
    }
 
    // Output:
    // println!("{:?}", representatives);
    println!("Generated {} clusters", total);
    println!("{} clusters contain more than one read", g1);
    println!("{} alignment calls", aln_calls);
    return (clusters, representatives, minimizer_database);
}

fn get_seq_hpol_comp<'a> (sequence: &'a [u8]) -> (Vec<u8>, Vec<usize>) {
    let mut last_val: u8 = 0; // TODO, check if safe.
    let mut current_count = 0;
    let mut all_read_hpol_lengths = Vec::new();
    let seq_hpol_comp: Vec<_> = sequence.iter().filter(
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

    return (seq_hpol_comp, all_read_hpol_lengths);
}

fn get_best_cluster_block_align(read_cl_id: usize, representatives: &HashMap<usize, InputSequence>, hit_clusters_ids: &HashMap<usize, u64>, hit_clusters_hit_positions: &HashMap<usize, Vec<usize>>, args_k: i64, args_aligned_threshold: f64) 
-> (i64, usize, f64, Option<Vec<u8>>, Option<Vec<u8>>, f64) {
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

        let match_id_tailored = ((1.0 - error_rate_sum) * args_k as f64) as i64;
        
        let (s1_alignment, s2_alignment, alignment_ratio) = parasail_block_alignment(seq, c_seq, args_k, match_id_tailored, gap_opening_penalty);
      
        if alignment_ratio > args_aligned_threshold {
            return (*cl_id as i64, nm_hits, error_rate_sum, Some(s1_alignment), Some(s2_alignment), alignment_ratio);
        }
    }
    return (-1, 0, 0.0, None, None, 0.0);
}


fn parasail_block_alignment(s1: &[u8], s2: &[u8], k: i64, match_id: i64, opening_penalty: i32) ->
    (Vec<u8>, Vec<u8>, f64) {

    let mismatch_penalty = -2;
    let gap_ext = 1;
    let match_score = 2;

    // def parasail_block_alignment(s1, s2, k, match_id, match_score = 2, mismatch_penalty = -2, opening_penalty = 5, gap_ext = 1):

    let user_matrix = Matrix::create("ACGT", match_score, mismatch_penalty);
    let result = parasailors::semi_global_alignment_trace_scan_sat_cigar(s1, s2, opening_penalty, gap_ext, &user_matrix);
    let (s1_alignment, s2_alignment) = cigar_to_seq(&result.cigar_trace, s1, s2);

    let match_vector: Vec<i64> = s1_alignment.iter().zip(s2_alignment.iter()).map(|(n1,n2)| if n1 == n2 { return 1; } else { return 0; }).collect();
    let mut match_window: VecDeque<i64> = match_vector[0..k as usize].iter().cloned().collect();
    let mut current_match_count: i64 = match_window.iter().sum();

    let mut aligned_region = Vec::new();
    if current_match_count >= match_id {
        aligned_region.push(1);
    } else {
        aligned_region.push(0);
    }

    for new_m_state in match_vector[k as usize..].iter() {
        let prev_m_state = match_window.pop_front().unwrap();
        current_match_count = current_match_count - prev_m_state + new_m_state;

        match_window.push_back(*new_m_state);

        if current_match_count >= match_id {
            aligned_region.push(1);
        } else {
            aligned_region.push(0);
        }
    }

    let alignment_ratio = aligned_region.iter().sum::<i64>() as f64 / s1.len() as f64;

    return (s1_alignment, s2_alignment, alignment_ratio);

}

fn cigar_to_seq(cigar: &String, query: &[u8], refs: &[u8]) -> (Vec<u8>, Vec<u8>) {

    let mut cigar_tuples = Vec::new();
    let re = Regex::new(r"[=DXSMI]+").unwrap();
    let result: Vec<&str> = re.split(cigar).collect();
    let mut i = 0;

    for index in 0..(result.len() - 1) {
        let length = result[index];
        i += length.len();
        let type_ = cigar.as_bytes()[i] as char;
        i += 1;
        cigar_tuples.push((length.parse::<usize>().unwrap(), type_));
    }

    let mut r_index = 0;
    let mut q_index = 0;
    let mut q_aln = Vec::new();
    let mut r_aln = Vec::new();

    for (length_, type_) in cigar_tuples {
        if type_ == '=' || type_ == 'X' {
            q_aln.extend_from_slice(&query[q_index..(q_index + length_)]);
            r_aln.extend_from_slice(&refs[r_index..(r_index + length_)]);

            r_index += length_;
            q_index += length_;
        } else if type_ == 'I' {
            // insertion w.r.t. reference
            r_aln.extend_from_slice(&std::iter::repeat("-").take(length_).collect::<String>().as_bytes());
            q_aln.extend_from_slice(&query[q_index..(q_index + length_)]);
          
            q_index += length_;

        } else if type_ == 'D' {
            // deletion w.r.t. reference
            r_aln.extend_from_slice(&refs[r_index..(r_index + length_)]);
            q_aln.extend_from_slice(&std::iter::repeat("-").take(length_).collect::<String>().as_bytes());
   
            r_index += length_;
        } else {
            println!("Error!");
            println!("{}", cigar);
            panic!();
        }
    }
    return (q_aln, r_aln);
}

fn get_best_cluster(read_cl_id: usize, compressed_seq_len: usize, hit_clusters_ids: &HashMap<usize, u64>, hit_clusters_hit_positions: &HashMap<usize, Vec<usize>>, minimizers: &Vec<(&[u8], usize)>, nummber_of_minimizers: usize, hit_clusters_hit_index: &HashMap<usize, Vec<usize>>, representatives: &HashMap<usize, InputSequence>, p_emp_probs: &HashMap<String, f64>, args_min_shared: usize, args_min_fraction: f64, args_min_prob_no_hits: f64, args_mapped_threshold: f64)
    -> (i64, usize, f64) {

    /*
          Tally up total covered (mapped region) and compare it with total unmapped region. What is counted as a consecutive mapped block depends on minimizer qualitity
        threshold as mapped-read-legnth/total-read-length is chosen to classify a read as belonging to the cluster
        
        Return: An integer >= 0 denoting the cluster ID that this read was assigned to. In not assigend to any previous cluster, return None. 
                [Also returing mapped ratio and nr shared minimizers to best read for logging purposes.]
    */

    let nr_shared_kmers = 0;
    let mut mapped_ratio = 0.0;

    if hit_clusters_ids.len() == 0 {
        return (-1, nr_shared_kmers, mapped_ratio);
    }

    let mut top_matches: Vec<(&usize, &Vec<usize>)> = hit_clusters_hit_positions.iter().collect();
    top_matches.sort_by_key(|x| (x.1.len(), x.1.iter().sum::<usize>(), representatives.get(x.0).unwrap().acc));
    top_matches.reverse();

    // println!("{:?}", top_matches);
    let top_hits = top_matches.get(0).unwrap().1.len();
    let nr_shared_kmers = top_hits;

    // println!("{}\t{}", top_hits, nr_shared_kmers);
    if top_hits < args_min_shared {
        // pass
    } else {
        for tm in top_matches {
            let cl_id = tm.0;
            let nm_hits = tm.1.len();
            // println!("{} <? {}", nm_hits, args_min_fraction * (top_hits as f64));
            if (nm_hits as f64) < args_min_fraction * (top_hits as f64) || nm_hits < args_min_shared {
                break;
            }

            let minimizer_hit_positions = hit_clusters_hit_positions.get(cl_id).unwrap();
            // println!("{:?}", minimizer_hit_positions);
            let minimizer_hit_indices = hit_clusters_hit_index.get(cl_id).unwrap();
            assert_eq!(minimizer_hit_indices.len(), minimizer_hit_positions.len(), "Unequal lengths!");

            let error_rate_c = representatives.get(cl_id).unwrap().error_rate;
            let error_rate_read = representatives.get(&read_cl_id).unwrap().error_rate;

            let p_error_in_kmers_emp = 1.0 - p_shared_minimizer_empirical(error_rate_read, error_rate_c, p_emp_probs);
            let minimizer_error_probabilities = vec![p_error_in_kmers_emp; nummber_of_minimizers];
            let mut total_mapped = 0;

            let mut prob_all_errors_since_last_hit = Vec::new();
            prob_all_errors_since_last_hit.push(minimizer_error_probabilities[..minimizer_hit_indices[0]].iter().fold(1.0, |a, b| a*b));
            for (hit_idx1, hit_idx2) in minimizer_hit_indices[..minimizer_hit_indices.len() - 1].iter().zip(minimizer_hit_indices[1..].iter()) {
                if *hit_idx1 == *hit_idx2 {
                    prob_all_errors_since_last_hit.push(1.);
                } else {
                    prob_all_errors_since_last_hit.push(minimizer_error_probabilities[*hit_idx1 + 1..*hit_idx2].iter().fold(1.0, |a, b| a*b));
                }
            }
            prob_all_errors_since_last_hit.push(minimizer_error_probabilities[minimizer_hit_indices[minimizer_hit_indices.len() - 1]+1..].iter().fold(1.0, |a, b| a*b));

            assert_eq!(prob_all_errors_since_last_hit.len(), minimizer_hit_positions.len() + 1, "Invalid lengths!");

            for i in 0..minimizer_hit_indices.len() {
                // println!("{} <? {}: {}", prob_all_errors_since_last_hit[i], args_min_prob_no_hits, prob_all_errors_since_last_hit[i] < args_min_prob_no_hits);

                if prob_all_errors_since_last_hit[i] < args_min_prob_no_hits {
                    // pass
                } else {
                    if i == 0 {
                        total_mapped += minimizer_hit_positions[i];
                    } else {
                        // println!("{}\t{}", minimizer_hit_positions[i], minimizer_hit_positions[i-1]);
                        total_mapped += minimizer_hit_positions[i] - minimizer_hit_positions[i - 1];
                    }
                }
            }
            // println!("{}", total_mapped);
            if prob_all_errors_since_last_hit[prob_all_errors_since_last_hit.len() - 1] < args_min_prob_no_hits {
                // pass
            } else {
                total_mapped += compressed_seq_len - minimizer_hit_positions[minimizer_hit_positions.len() - 1];
            }

            mapped_ratio = total_mapped as f64 / compressed_seq_len as f64;

            if mapped_ratio > args_mapped_threshold {
                return (*cl_id as i64, nm_hits, mapped_ratio);
            }
            // [reduce(mul, minimizer_error_probabilities[: minimizer_hit_indices[0]], 1)] +  
            // [ reduce(mul, minimizer_error_probabilities[hit_idx1+1: hit_idx2], 1) for hit_idx1, hit_idx2 in zip(minimizer_hit_indices[:-1], minimizer_hit_indices[1:]) ] + 
            // [reduce(mul, minimizer_error_probabilities[minimizer_hit_indices[-1]+1 : ], 1)]

        }
    }
    
    return (-1, nr_shared_kmers, mapped_ratio);
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

fn get_all_hits(minimizers: &Vec<(&[u8], usize)>, minimizer_database: &HashMap<Vec<u8>, Vec<usize>>, read_cl_id: usize) 
    -> (HashMap<usize, u64>, HashMap<usize, Vec<usize>>,HashMap<usize, Vec<usize>>) {

    //     Get all representatives ID's that shares matches with the minimizers in the read.  

    let mut hit_clusters_ids = HashMap::new();
    let mut hit_clusters_hit_index = HashMap::new();
    let mut hit_clusters_hit_positions = HashMap::new();

    // println!("{:?}", minimizers);
    for (i, (m, pos)) in minimizers.iter().enumerate() {
        if minimizer_database.contains_key(*m) {
            for cl_id in minimizer_database.get(*m).unwrap() {
                // println!("m: {:?}\tcl_id: {}", *m, cl_id);
                hit_clusters_ids.insert(*cl_id, hit_clusters_ids.get(&cl_id).unwrap_or(&0) + 1);
                hit_clusters_hit_index.entry(*cl_id).or_insert(Vec::new()).push(i);
                hit_clusters_hit_positions.entry(*cl_id).or_insert(Vec::new()).push(*pos);
                // println!("{:?}", pos)
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
    for i in 0..(w + 1) as usize {
        if i + (k_size as usize) < seq.len() {
            window_kmers.push_back(&seq[(i as usize)..((i + k_size as usize) as usize)]);
        }
    }

    let mut minimizers = Vec::new();
    let mut curr_min = get_min_and_index(&window_kmers);
    minimizers.push(curr_min);

    // println!("{}", String::from_utf8_lossy(curr_min.0));
    for i in (w as usize+1)..(seq.len() - k_size as usize + 1) {
        let new_kmer = &seq[i..(i + k_size as usize)];
        let discarded = window_kmers.pop_front().unwrap();
        window_kmers.push_back(new_kmer);

        if curr_min.0 == discarded {
            curr_min = get_min_and_index(&window_kmers);
            curr_min.1 += i - w as usize;
            minimizers.push(curr_min);
            // println!("Pushing {:?}", curr_min);
        } else if new_kmer < curr_min.0 {
            curr_min = (new_kmer, i);
            minimizers.push(curr_min);
            // println!("Pushing {:?}", curr_min);
        }
    }
    // println!("{:?}", minimizers);
    return minimizers;
}


