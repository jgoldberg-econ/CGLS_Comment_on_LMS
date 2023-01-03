# README for replication of "A Comment on: Low Interest Rates, Market Power, and Productivity Growth," by Craig A. Chikis, Jonathan Goldberg, and David López-Salido (2022).


## How to run code 
Clone the Git repo. Set your working directory to the repo and run `main.m`.


## Approach

Our code heavily leverages the LMS replication code, downloaded from the Econometrica website, so that a reader can see clearly and precisely how we deviate from their code.  Comments (%) in the scripts and functions indicate where code is obtained directly from LMS's code without modification.  Where the LMS code is modified, comments explain the modifications.  

The MATLAB script `transition_figs.m` replicates Figures 1 and 2 in our comment and Figure IA.2 of our Online Appendix.

The MATLAB script `bgp_figs.m` replicates Figure 3 in our comment and Figure IA.3 of our Online Appendix.  (Figure IA.1 of the Online Appendix is obtained using the methods described in our Online Appendix A.5.)

## Previous draft

For a previous version of our comment, two complete and distinct replication codes are provided. See https://github.com/cchikis/CommentOnLMS_CGLS.  One of the replication codes uses "standard" numerical methods from the Aghion-Howitt literature.  Specifically, we solve for investment success rates using uniformization and value function iteration (see Acemoglu and Akcigit, JEEA 2012).  This code allows a reader to see the HJB equations in an intuitive way.  This code also stores all rates not in percent. 
