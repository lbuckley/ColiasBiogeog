LambdaModel output: 
Copoints.csv: grid cell information

lambda1_....abs: fitness and related predictions (see below)
dimensions:1: years2: grid cells3: absorptivity4: generation4 files contain output for Lambda, FAT (flight activity times, hours), Egg Viability (proportion), body temperature (C)

PupTemps_...rds: pupal temperature and related estimates 
dimensions:
1: metrics: "stat","yr","gen","Jlarv", "Jpup","Jadult","Tlarv","Tpup","Tad","Tlarv_fixed","Tpup_fixed","Tad_fixed"2: years3: grid cells4: generation 

absmean.abs: absorptivity predictions from evolutionary model
dimensions:
1: years 
2: grid cells
3: generations
4: evolutionary scenarios :no plast, plast, only plast
5: metrics: abssample, absmid, rn, Babsmid, Brn)

Lambda mean.abs: fitness predictions from evolutionary model
dimensions:
1: years 
2: grid cells
3: generations
4: evolutionary scenarios :no plast, plast, only plast

Abs.opt_...rds: optimal absorptivity, 
dimensions: years, grid cells, generation

See ColiasBiogeog_EvolModel_June2017_newFigures.R for formatting information and to make plots.
