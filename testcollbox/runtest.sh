# small-scale test of collision module
# Particles in a periodic box are allowed
# to bounce elastically. Start with uniform
# speeds, will eventually evolve into a
# Maxwell-Bolztmann distribution.

../charmrun +p 2 ../ChaNGa -v 1 testcollbox.param ++local
