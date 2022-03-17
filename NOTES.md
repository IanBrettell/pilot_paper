# Notebook

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

