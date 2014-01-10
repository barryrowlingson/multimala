\name{read_multimala_output}
\alias{read_multimala_output}
\title{Read Multimala outputs}
\usage{
  read_multimala_output(filename, grid, mask)
}
\arguments{
  \item{filename}{one of the output files}

  \item{grid}{grid specification}

  \item{mask}{optional mask}
}
\value{
  a raster of the values from the file
}
\description{
  Read multimala output files
}
\details{
  multimala generates simple streamed output files. This
  loads them into raster objects
}
\author{
  Barry S Rowlingson
}

