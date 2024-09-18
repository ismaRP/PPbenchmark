library(AmyloGram)
library(appnn)
blg = read_txt('~/palaeoproteomics/BLG/bovin_blg/bovin_blgA.fasta')
amylogram_pred = predict(AmyloGram_model, blg)
amylogram_pred = amylogram_pred$detailed
amylogram_pred$position = 1:178
write.csv(amylogram_pred,
          '~/palaeoproteomics/BLG/bovin_blg/amyloid_propensity/AmyloGram.csv',
          row.names = F, quote = F)


appnn_pred = appnn(paste(blg[[1]], collapse = ''))
appnn_pred = data.frame(
  position = 1:178,
  aminoacid = blg[[1]],
  prediction = appnn_pred[[1]]$aminoacids
)
write.csv(appnn_pred,
          '~/palaeoproteomics/BLG/bovin_blg/amyloid_propensity/appnn_amylopred.csv',
          row.names = F, quote = F)
