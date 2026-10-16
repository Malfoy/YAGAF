pub mod consensus;
pub mod index;
pub mod input;
pub mod pipeline;

#[derive(Debug)]
pub struct Scaffold {
    pub header: Vec<u8>,
    pub sequence: Vec<u8>,
}

#[derive(Debug)]
pub struct Gap {
    pub scaffold: usize,
    pub start: usize,
    pub end: usize,
    pub anchors: Option<[Vec<u8>; 2]>,
    pub rejection: Option<&'static str>,
}

pub fn is_base(base: u8) -> bool {
    matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

pub fn revcomp(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

pub fn discover(scaffolds: &[Scaffold], k: usize) -> Vec<Gap> {
    let mut gaps = Vec::new();
    for (scaffold, record) in scaffolds.iter().enumerate() {
        let seq = &record.sequence;
        let mut pos = 0;
        while pos < seq.len() {
            if !seq[pos].eq_ignore_ascii_case(&b'N') {
                pos += 1;
                continue;
            }
            let start = pos;
            while pos < seq.len() && seq[pos].eq_ignore_ascii_case(&b'N') {
                pos += 1;
            }
            let end = pos;
            let anchors = start.checked_sub(k).and_then(|left| {
                let right_end = end.checked_add(k)?;
                if right_end > seq.len() {
                    return None;
                }
                let a = &seq[left..start];
                let b = &seq[end..right_end];
                (a.iter().chain(b).all(|&base| is_base(base)))
                    .then(|| [a.to_ascii_uppercase(), b.to_ascii_uppercase()])
            });
            let rejection = anchors.is_none().then_some("missing_or_ambiguous_flank");
            gaps.push(Gap {
                scaffold,
                start,
                end,
                anchors,
                rejection,
            });
        }
    }
    gaps
}
