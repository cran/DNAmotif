\name{DNAmotifs}
\alias{DNAmotifs}
\title{Generation of motifs from DNA sequences}
\usage{DNAmotifs(fasta_file, ws, cut_off)}
\description{
Using a fasta file as input, the motifs and its corresponding frequencies are generated, considering threshold values (window size and similarity), by making consensus segments via local alignment with gap.
}
\arguments{
\item{fasta_file}{Sequence file path (.fasta format)}
\item{ws}{Window size}
\item{cut_off}{Minimum similarity percentage between the motifs for generating a consensus motif}
}
\value{
\code{final_results} A dataframe of motifs and their corresponding frequncies
}
\examples{
example_fasta = system.file("exdata/sample.fasta", package = "DNAmotif")
DNAmotifs(fasta_file = example_fasta, ws = 15, cut_off = 0.75)
}
