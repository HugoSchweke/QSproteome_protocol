## arg must be model id + suffix _all_csym.dat, ie UNIPROT_VX_Y_all_csym.dat
args = commandArgs(trailingOnly=TRUE)


#path = "/data5/elevy/01_3dcomplexV0/results/005_syms/ananas_clashscore_allsymC_nodiso80/"
file = args[1]
print(file)

thr_clashes = 200

list_pdb = list()
#for (file in clashfiles) {
code = gsub("_all_csym.dat", "", basename(file))
#dat = read.table(paste0(path,file), header = T, fill = T, stringsAsFactors = F)
dat = read.table(paste0(file), header = T, fill = T, stringsAsFactors = F)
dat = na.omit(dat)
dat = dat[dat$clashscore<thr_clashes,] # Remove all clash scores above 100, considered as faulty structures
dat = dat[order(dat$av.rmsd),] # Rank by rmsd, we just take best rmsd among the structures filtered on clashscore
list_pdb[[code]] = dat
##print(list_pdb[[code]])
#}

best.rmsd = lapply(list_pdb, 
                   function(x) {
                      if (nrow(x) >= 1) {
                        bestrmsd = x[1,2]
                      } else {
                        bestrmsd = NA
                      }
                      return(bestrmsd)
                   })
best.rmsd = unlist(best.rmsd)

best.sym = lapply(list_pdb, 
                  function(x) {
                      print(x)
                    if (nrow(x) >= 1) {
                      best.sym = x[1,1]
                    } else {
                      best.sym = "NPS"
                    }
                    return(best.sym)
                  })
best.sym = unlist(best.sym)
print(best.sym)

best.clash = lapply(list_pdb, 
                    function(x) {
                      if (nrow(x) >= 1) {
                        best.clash = x[1,3]
                      } else {
                        best.clash = NA
                      }
                      return(best.clash)
                    })
best.clash = unlist(best.clash)

data.best = data.frame(code = names(best.clash),
                       symmetry = best.sym,
                       rmsd = best.rmsd,
                       clash.score = best.clash)
print(data.best)

write.csv(data.best, file = paste0("/data5/elevy/01_3dcomplexV0/results/005_syms/ananas_clashscore_bestsym_nodiso80_thr400/", data.best$code, "_best_sym_clash.csv"), quote = F, row.names = F)
