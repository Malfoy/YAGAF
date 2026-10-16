use crate::graphs::poa::POAGraph;
use petgraph::{dot::Dot, graph::IndexType};
use std::fmt;
pub fn format_as_dot<Ix: IndexType>(
    writer: &mut impl fmt::Write,
    graph: &POAGraph<Ix>,
) -> fmt::Result {
    let transformed = graph.graph.map(
        |ix, data| format!("{:?} ({:?})", char::from(data.symbol), ix.index()),
        |_, data| format!("{}, {:?}", data.weight, data.sequence_ids),
    );

    let dot = Dot::new(&transformed);

    writeln!(writer, "{}", dot)?;

    Ok(())
}
