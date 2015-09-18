# Author: Babak Naimi, naimi.b@gmail.com
# Date :  August 2014
# Version 1.0
# Licence GPL v3 

setClassUnion('data.frameORmatrix',c("data.frame","matrix"))
setClassUnion('matrixORnull',c("NULL","matrix"))

setClass("Entrogram",
         representation(width="numeric",
                        cutoff="numeric",
                        entrogramCloud='matrixORnull',
                        entrogram="data.frame"
                        ),
         prototype(
           entrogramCloud=NULL
           )
)
