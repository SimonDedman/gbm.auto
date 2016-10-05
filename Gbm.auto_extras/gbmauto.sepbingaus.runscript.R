gbm.auto2(samples = Juveniles,
          grids = mygrids,
          expvar = c(4:10,14,16,20,24,28,36),
          resvar = 43,
          tc = c(2,13),
          lr = c(0.01,0.005),
          bf = 0.5,
          mapshape = Crop_Map)

gbm.auto2(samples = Juveniles,
          grids = mygrids,
          expvar = c(4:10,14,16,20,24,28,36),
          resvar = 43,
          tc = list(c(2,13),2),
          lr = list(c(0.01,0.005),0.0001),
          bf = list(c(0.5,0.7),0.5),
          mapshape = Crop_Map)

gbm.auto2(samples = Juveniles,
          grids = mygrids,
          expvar = c(4:10,14,16,20,24,28,36),
          resvar = 43,
          tc = list(2, c(2,13)),
          lr = list(0.01, c(0.01,0.005)),
          bf = 0.5,
          mapshape = Crop_Map, simp = FALSE, varint = FALSE, map = FALSE, BnW = FALSE, RSB = F, linesfiles = FALSE)

gbm.auto2(samples = Juveniles,
          +           grids = mygrids,
          +           expvar = c(4:10,14,16,20,24,28,36),
          +           resvar = 43,
          +           tc = c(2,3),
          +           lr = c(0.01,0.005),
          +           bf = c(0.5,0.7),
          +           mapshape = Crop_Map, simp = FALSE, varint = FALSE, map = FALSE, BnW = FALSE, RSB = F, linesfiles = FALSE, savegbm = F)
[1] 2 3
[1] 2 3
[1] 0.010 0.005
[1] 0.010 0.005
[1] 0.5 0.7
[1] 0.5 0.7

gbm.auto2(samples = Juveniles,
          grids = mygrids,
          expvar = c(4:10,14,16,20,24,28,36),
          resvar = 43,
          tc = 2,
          lr = list(c(0.01,0.02),0.001),
          bf = 0.5,
          mapshape = Crop_Map, simp = FALSE, varint = FALSE, map = FALSE, BnW = FALSE, RSB = F, linesfiles = FALSE, savegbm = F)
