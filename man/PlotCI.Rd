\name{PlotCI}
\alias{PlotCI}
\title{
Plot Simultaneous Confidence Intervals
}
\description{
Draw a plot of simultaneous confidence intervals (SCIs) from multiple comparisons in a longitudinal setting using the \code{ggplot2} framework.
}
\usage{
PlotCI(output, title="95\% SCIs")
}
\arguments{
  \item{output}{An object of class \code{silo} obtained from using functions \code{SimLongi} or \code{SimLongiMMM}, or a named list of (up to seven) such objects.}
  \item{title}{An optional character string giving the plot title. By default, a standard headline ("95\% SCIs") is printed.}
}
\details{
If \code{output} is a list containing several objects of class \code{silo}, this list must be named (e.g., \code{list(name1=firstobject, name2=secondobject)}) so that a legend can be created according to the names.

If no headline is desired, simply enter \code{title=NULL}.
}
\value{
A graphical display of SCIs.
}
\author{
Philip Pallmann \email{pallmann@biostat.uni-hannover.de}
}
\seealso{
\code{\link{SimLongi}}, \code{\link{SimLongiMMM}}
}
\examples{
data(heart)

# Many-to-one comparisons of groups per time point
# taking the third group ("control") as reference

mycomp <- SimLongi(data=heart, response="heartrate", group="drug",
                   time="time", id="person", direction="gpt",
                   type="Dunnett", base=3, df="adj")
      
PlotCI(mycomp)

mycompMMM <- SimLongiMMM(data=heart, response="heartrate", group="drug",
                         time="time", id="person", type="Dunnett", base=3,
                         df="adj")
                   
PlotCI(list(SL=mycomp, SL_MMM=mycompMMM))
}
\keyword{hplot}