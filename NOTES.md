# Notebook

## 29 April 2022

Ewan's ideas for figures:

* Figure 1 diagram of the set up, the pairing and then some tracks

* Figure 2 is something about “our tracking is good”

* Figure 3 polar plots and boxplots for DGE, plus HMM tiles

* Figure 4 is the same for SGE

* One thought on the tracks, there is a big time dependence of some of the states (most obviously the still states at the start) - I was wondering if we should pick a couple of states which show strongest time dependence (and / or time by strain dependence)

## 16 March 2022

Tracking success:

* 188/614 (30%) perfectly tracked

* 483/614 (78%) > 99% tracked

* Minimum is 85%

Progress to data analysis.

## 20220311

Tracking success: 

* 120/614 (20%) perfectly tracked

* 383/614 (62%) > 99% tracked

Most tracking errors caused by lack of background subtraction.
Make bgsub = True as the default, as these inbred lines tend to move around a lot relative to MIKK, which makes the "medaka ghost" phenomenon when bgsub = True less likely.

