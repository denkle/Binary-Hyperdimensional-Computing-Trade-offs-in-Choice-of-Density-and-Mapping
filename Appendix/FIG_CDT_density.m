%This script outputs the figure corresponding to Fig. 12 in the article

DENS=[0.01]; % density of vectors
computeCDT_experimental_theoretical(DENS)
DENS=[0.02]; % density of vectors
computeCDT_experimental_theoretical(DENS)
DENS=[0.03]; % density of vectors
computeCDT_experimental_theoretical(DENS)

text(2.5,0.35,'p_1=0.03')
text(2.5,0.2,'p_1=0.02')
text(2.5,0.07,'p_1=0.01')

xlabel('Number of iterations in CDT, \it{T}')
ylabel('Density of ones, p_{T1}')
box on
