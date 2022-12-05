# PI3K-PTEN-model
A program to generate the figures in the paper  http://doi.org/10.1242/jcs.108373

Shibata, T., Nishikawa, M., Matsuoka, S., & Ueda, M. (2012). Modeling the self-organized phosphatidylinositol lipid signaling system in chemotactic cells using quantitative image analysis. Journal of Cell Science, 125(21), 5138–5150. http://doi.org/10.1242/jcs.108373

### Files
- PI3KPTENModel.c is a c-code program, which performs the numerical calculation of Eq. (4) in the paper.
- param_A.txt to param_E.txt give the parameters in Fig. 7 for the program.
- makeKymograph.m is a matlab program to generate a kymograph a data file.

### How to run
- `# make`
- `# PI3KPTENModel.out -f param_A.txt -o outA`
- `# for i in {A..E};do ./PI3KPTENModel.out -o out${i} -f param_${i}.txt done`
- `# matlab -nodesktop -nosplash << 'EOF`
- `makeKymograph;`
- `exit`
- `EOF`
