## Introduction


### Communication from B.B. to E.M.


B.B. to do CellHashing.R and SeuratFunctions.R with the goal of: 

- correctly tagging the minimal set for @export, 

- adding @import/@importFrom to get it working.  

- minimal comment block for anything that is not exported, if anything.  We want some documentation, but in order to reduce work, focus on the primary exported ones.  


3)  import() can be used instead of importFrom(); however, it seems like we ultimately want to be using importFrom to be more surgical.  

  we need to be using import() or importFrom()..comments? 
  
  https://kbroman.org/pkg_primer/pages/docs.html
  
4) When I run devtools::test() it keeps failing b/c it cant find various methods that should be loaded.

5)  convert function names to capitalized

6) save the cell cycle gene lists in the package

7) get  install_github for OOSAP working


