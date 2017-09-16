# Abstract
<div style="text-align: justify; text-indent: 30px;"><p>
Recently Magnetic Resonance Imaging (MRI) Compressive Sensing (CS) reconstruction algorithms attempt to incorporate adaptive dictionaries to enhance reconstruction performance. However, a fixed-dictionary baseline is missing which prevents a fair evaluation of the role of adaptive dictionaries in reconstruction. Herein, a Fixed-Dictionary Baseline (FDB) algorithm is proposed, which is a fast, easy-to-use, and strong algorithm based on the latest MRI CS algorithms. It is compared with the current baseline, DLMRI, and state-of-the-art algorithms, UNITE-MRI and DDTF-MRI under various conditions including dictionary redundancy. Experimentally, the redundancy considerably enhanced all algorithms (up to 7dB PSNR with FDB). Under radial down sampling, DDTF-MRI achieves 1 dB PSNR gain over FDB only at up to 3-fold acceleration; UNITE-MRI and DLMRI show no gain. Under 1D down sampling, FDB outperforms all the adaptive-dictionary algorithms up to 0.9 dB PSNR. Overall, DLMRI does considerably worse than the other three, including FDB. Our study reveals that, the previous baseline DLMRI is inferior to the simpler FDB algorithm. Reporting marked improvement over DLMRI can thus be misleading. Adaptive dictionaries can enhance reconstruction when acceleration is low, but their ability to push the speed limit remains to be proven. With high redundancy, the benefit of adaptive dictionary tends to diminish.
</p></div>

# Results
<center>
<img width="440" height="1125" align="center" src="assets/images/reconstruction_comparison.png?raw=true">
</center>
<div style="text-align: justify; text-align:center;" >
<p>
Reconstruction with radial measurements with 5-fold acceleration.
(a) is the original image. (b) is the pseudo radial down sampling pattern.
(c)(e)(g)(i) show the reconstructed images using FDB, DDTF-MRI, UNITE-MRI,
and DLMRI, with maximum block/generator overlapping.
(d)(f)(h)(j) show the corresponding error maps. Reconstruction
PSNR are 39.1, 39.2, 38.0, and 32.4 respectively.
The enhancement of DDTF-MRI and UNITE-MRI over DLMRI agree with previous reports.
</p></div>

# Publications
<div style="text-align: justify;"><p>
Yuan Wang, Yao Wang, and Yvonne W Lui,
A new fixed-dictionary baseline for compressive sensing magnetic resonance imaging reconstruction algorithms using adaptive dictionaries,
submitted to IJBI, under review.
</p></div>

# code usage
<div style="text-align: justify;"><p>
Start with and walk through the demo.m,
and you should be able to get an idea how it works.
</p></div>
> At the beginning the demo.m, you need to set your Matlab working directory to the 'code'.

#### <a href="https://yuanwangonward.github.io/">Back to Yuan's homepage</a>

