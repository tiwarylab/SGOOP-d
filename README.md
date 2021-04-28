# SGOOP-d: Estimating reaction coordinate dimensionality and kinetic distances for rare event systems


Here we show How we do a single step of SGOOP. Note that the off-diagonal terms are strictly positive, so the diagonal terms are all negative. Therefore, the eigenvalues of the rate matrix we define here are negative. The transition rate associated to each eigenmode is **negative eigenvalues**.

<img src="https://user-images.githubusercontent.com/22850008/115948672-6ec29a00-a49d-11eb-81e2-1935b6d40ff9.png" width="550">

SGOOP-d further utilizes the eigenvalues (rate) and eigenvectors (RC) from the MaxCal-derived rate matrix of each multi-SGOOP RC. To estimate the timescale in each eigenspace, we use absolute value of the inverse eigenvalues.

<img src="https://user-images.githubusercontent.com/22850008/116347501-5236a280-a7ba-11eb-8254-f8cbcb61c636.png" width="700">
