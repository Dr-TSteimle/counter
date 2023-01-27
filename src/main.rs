use rust_htslib::bam::{IndexedReader, Read};
use std::fs::File;
use std::io::BufReader;
use csv::ReaderBuilder;
use std::time::Instant;
use clap::Parser;

// find /Turbine-pool/LAL-T_RNAseq/BAM_star/*.sorted.bam | parallel -j30 /home/thomas/Documents/Programmes/counter/target/debug/counter {} NGS/tools/arriba_v2.1.0/RefSeq_hg38.gtf ">" {}.count
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
   #[arg(short='a', long)]
   bam_path: String,

   #[arg(short, long)]
   bed_path: String,
}

struct Position {
    contig: String,
    from  : i32,
    to    : i32
}
fn main() {
    let args = Args::parse();
    let path_bam = args.bam_path; // "/Turbine-pool/LAL-T_RNAseq/BAM_star/58_MAS.star_out.sorted.bam"
    let path_bed = args.bed_path; // "/home/thomas/NGS/tools/arriba_v2.1.0/RefSeq_hg38.gtf"
    let mut bam = IndexedReader::from_path(&path_bam).unwrap();
    count_on(&path_bed, &mut bam);
    std::process::exit(0);
}

fn get_from_pos(pos: Position, bam: &mut IndexedReader) -> i128 {
    eprintln!("Reading position {}:{}-{}", pos.contig, pos.from, pos.to);
    let mut n_reads: i128 = 0;
    let r = bam.fetch((&pos.contig, pos.from, pos.to));
    
    for read in bam.records() {
        if read.is_ok() {n_reads += 1;}
    }
    n_reads
}

fn count_on(path: &String, bam: &mut IndexedReader) {
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

            let res = get_from_pos(Position {contig: line[0].clone(), from: line[1].parse().unwrap(), to: line[2].parse().unwrap()}, bam);
            line.push(res.to_string());
            println!("{}", line.join(&"\t"));           

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