use rust_htslib::bam::{IndexedReader, Read};
use std::fs::File;
use std::io::BufReader;
use csv::ReaderBuilder;
use std::time::Instant;
use regex::Regex;

// find /Turbine-pool/LAL-T_RNAseq/BAM_star/*.sorted.bam | parallel -j30 /home/thomas/Documents/Programmes/counter/target/debug/counter {} NGS/tools/arriba_v2.1.0/RefSeq_hg38.gtf ">" {}.count

struct Position {
    contig: String,
    from: i32,
    to: i32
}
fn main() {
    let path_bam = std::env::args().nth(1).expect("no bam path given"); // "/Turbine-pool/LAL-T_RNAseq/BAM_star/58_MAS.star_out.sorted.bam"
    let path_gtf = std::env::args().nth(2).expect("no gtf path given"); // "/home/thomas/NGS/tools/arriba_v2.1.0/RefSeq_hg38.gtf"
    let mut bam = IndexedReader::from_path(&path_bam).unwrap();
    count_on(&path_gtf, &mut bam, Regex::new(r"RNA[0-9]{1,2}|RNU[0-9]|RNVU[1-9]").unwrap());
}

fn get_from_pos(pos: Position, bam: &mut IndexedReader) -> (i32,i32) {
    let mut lane_counts: (i32,i32) = (0,0);
    let _ = bam.fetch((&pos.contig, pos.from, pos.to));

    for read in bam.records() {
        let mut n_sep: i32 = 0;
        let lane1: u8 = '1' as u8;
        let lane2: u8 = '2' as u8;
        let read = read.unwrap();
        read.qname().to_vec().iter().for_each(|a| {
            if n_sep == 3 {
                if *a == lane1 {
                    lane_counts.0 += 1;
                } else if *a == lane2 {
                    lane_counts.1 += 1;
                }
            }
            if *a == b":"[0] {
                n_sep+=1
            }
        });
    }
    lane_counts
}

fn count_on(path: &String, bam: &mut IndexedReader, omit_line: Regex) {
    let mut n_iter = 1;
    let mut n_reads = 1;
    let now = Instant::now();

    let file: File = File::open(path).unwrap();
    let mut reader = ReaderBuilder::new()
        // .comment(Some(b'!'))
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(BufReader::new(file));
    
    for line in reader.records() {
        let mut line: Vec<String> = line.unwrap().iter().map(|a| a.to_string()).collect();
        if !omit_line.is_match(&line[8]) {
            let res = get_from_pos(Position {contig: line[0].clone(), from: line[3].parse().unwrap(), to: line[4].parse().unwrap()}, bam);
            line.push(res.0.to_string());
            n_reads += res.0;
            line.push(res.1.to_string());
            n_reads += res.1;
            println!("{}", line.join(&"\t"));           
        }
        n_iter += 1;
        if n_iter % 10_000 == 0 {
            eprintln!("{} positions parsed in {} seconds", n_iter, now.elapsed().as_secs());
        }
        if n_reads % 100_000 == 0 {
            eprintln!("{} reads parsed in {} seconds", n_reads, now.elapsed().as_secs());
        }
    }
    eprintln!("{} reads assigned to {} positions in {} seconds.", n_reads, n_iter, now.elapsed().as_secs());
}