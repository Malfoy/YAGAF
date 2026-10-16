use petgraph::graph::IndexType;
use petgraph::visit::EdgeRef;
use rustc_hash::FxHashMap as HashMap;
use rustc_hash::FxHashSet;
use serde::de::DeserializeOwned;
use std::cell::RefCell;
use std::io::Write;
use std::iter;

use crate::errors::PoastaError;
use crate::graphs::poa::{POAGraph, POANodeIndex};
use crate::graphs::AlignableRefGraph;

fn fasta_aln_for_seq<Ix>(
    graph: &POAGraph<Ix>,
    node_to_column: &HashMap<POANodeIndex<Ix>, usize>,
    seq_id: usize,
    start_node: POANodeIndex<Ix>,
    width: usize,
) -> Vec<u8>
where
    Ix: IndexType + DeserializeOwned,
{
    let mut seq = Vec::with_capacity(width);
    let mut curr = Some(start_node);
    // Columns are zero-based. Track the next output column, rather than using
    // zero as both a sentinel and the first aligned column (upstream bug).
    let mut next_col = 0;
    while let Some(n) = curr {
        let node_col = node_to_column[&n];

        let gap_length = node_col - next_col;
        seq.extend(iter::repeat(b'-').take(gap_length));
        seq.push(graph.get_symbol(n));

        curr = None;
        for out_edge in graph.graph.edges(n) {
            if out_edge
                .weight()
                .sequence_ids
                .binary_search(&seq_id)
                .is_ok()
            {
                curr = Some(out_edge.target())
            }
        }

        next_col = node_col + 1;
    }

    if next_col < width {
        let gap_length = width - next_col;
        seq.extend(iter::repeat(b'-').take(gap_length));
    }

    seq
}

pub fn poa_graph_to_msa<Ix>(graph: &POAGraph<Ix>) -> Result<Vec<Vec<u8>>, PoastaError>
where
    Ix: IndexType + DeserializeOwned,
{
    // First assign a column in the final MSA output for every node in the graph,
    // by visiting nodes in topological order, taking into account 'aligned_nodes'
    let mut node_to_column = HashMap::default();
    let mut stack = vec![(
        graph.start_node(),
        RefCell::new(Vec::from_iter(graph.successors(graph.start_node()))),
    )];

    let next_valid_child = |succ_iter: &RefCell<Vec<POANodeIndex<Ix>>>,
                            visited: &FxHashSet<POANodeIndex<Ix>>|
     -> Option<POANodeIndex<Ix>> {
        while let Some(succ) = succ_iter.borrow_mut().pop() {
            if !visited.contains(&succ) {
                return Some(succ);
            }
        }

        None
    };

    // DFS, visiting nodes in post order
    let mut visited = FxHashSet::default();
    let mut rev_postorder = vec![];
    while !stack.is_empty() {
        let (_, succ_iter) = stack.last().unwrap();

        if let Some(child) = next_valid_child(succ_iter, &visited) {
            visited.insert(child);
            let mut successors = Vec::from_iter(graph.successors(child));

            // When visiting a node, we also consider all aligned nodes visited,
            // since they will be part of the same column in the FASTA.
            // Let's also add their successors to the stack.
            for aln_node in graph.get_aligned_nodes(child) {
                if !visited.contains(aln_node) {
                    visited.insert(*aln_node);
                    successors.extend(graph.successors(*aln_node));
                }
            }
            stack.push((child, RefCell::new(successors)));
        } else {
            let (last, _) = stack.pop().unwrap();
            rev_postorder.push(last);
        }
    }

    rev_postorder.reverse();

    let mut curr_col = 0;
    for n in &rev_postorder {
        if *n == graph.start_node() || *n == graph.end_node() {
            continue;
        }

        if !node_to_column.contains_key(n) {
            node_to_column.insert(*n, curr_col);

            for aligned_node in graph.get_aligned_nodes(*n) {
                node_to_column.insert(*aligned_node, curr_col);
            }

            curr_col += 1;
        }
    }

    Ok(graph
        .sequences
        .iter()
        .enumerate()
        .map(|(seq_id, seq)| {
            fasta_aln_for_seq(graph, &node_to_column, seq_id, seq.start_node(), curr_col)
        })
        .collect())
}

pub fn poa_graph_to_fasta<Ix, W>(graph: &POAGraph<Ix>, mut output: W) -> Result<(), PoastaError>
where
    Ix: IndexType + DeserializeOwned,
    W: Write,
{
    for (sequence, row) in graph.sequences.iter().zip(poa_graph_to_msa(graph)?) {
        writeln!(output, ">{}", sequence.name())?;
        output.write_all(&row)?;
        output.write_all(b"\n")?;
    }
    Ok(())
}
