library(Seurat)

a = Sys.glob('/data/langenau/human_rms_pdxs/data_april/*Robj')

load(a[1])

system('mkdir -p ../results/normal_muscle')

seu = UpdateSeuratObject(AP10_SkM.Regr.S_G2M.Stress)
seu = as.loom(seu, filename=paste0('../results/normal_muscle/AP10_SkM.Regr.S_G2M.Stress.loom'), verbose=T)
seu$close_all()

load(a[2])

seu = UpdateSeuratObject(AP11.34_SkM.Regr.S_G2M.Stress)
seu = as.loom(seu, filename=paste0('../results/normal_muscle/AP11.34_SkM.Regr.S_G2M.Stress.loom'), verbose=T)
seu$close_all()

## head(AP10_SkM.Regr.S_G2M.Stress@meta.data)
## head(AP10_SkM.Regr.S_G2M.Stress@meta.data@Phase)
