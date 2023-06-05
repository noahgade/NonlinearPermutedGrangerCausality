import numpy as np
import rdata
import cmlp

dat = rdata.conversion.convert(rdata.parser.parse_file("simdata1000.RData"))

output = np.zeros((1600, 7, 7), dtype = float)

for dat_indx in range(1600):
    data = np.asmatrix(dat['simdata1000'][dat_indx])
    dimension = data.shape[-1]
    output[dat_indx, :, :] = cmlp.fit_cmlpGSGL(data, dimension, 3, 100, 0.5, 10000)
    
np.save("/home/ndg5e/NPGC2/cmlpGSGL1000.npy", output)