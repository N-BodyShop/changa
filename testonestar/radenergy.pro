function rad,m

if m>20 then return, 81.0*m^2.14*0.3029*m^(-2.7) $
else if m>2 then return, 1.78*m^3.5*0.3029*m^(-2.7) $
else return, 0.75*m^4.8*0.3029*m^(-2.7)
end

function kroupa,m
return, 0.3029*m^(-1.7)
end

print,qsimp('kroupa',1.0,40.0)
print,qsimp('rad',1.0,40.0)
end
