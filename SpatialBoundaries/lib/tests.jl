using Test

# Straight line

a = [1 1; 0 0]
@test _rateofchange(a) == (1.0, 0.0)

a = [1 0; 1 0]
@test _rateofchange(a) == (1.0, 90.0)

a = [0 0; 1 1]
@test _rateofchange(a) == (1.0, 180.0)

a = [0 1; 0 1]
@test _rateofchange(a) == (1.0, 270.0)

# Diagonal

a = [1 1; 0 1]
@test _rateofchange(a)[2] == 225

a = [1 0; 1 1]
@test _rateofchange(a)[2] == 45

a = [0 1; 1 1]
@test _rateofchange(a)[2] == 315

a = [1 1; 1 0]
@test _rateofchange(a)[2] == 135
