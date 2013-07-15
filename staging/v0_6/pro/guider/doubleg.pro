function doubleg, xx, sig1,sig2, amp1,amp2,comp1=comp1,comp2=comp2
   comp1 = amp1*exp(-xx*xx/(2*sig1*sig1))
   comp2 = amp2*exp(-xx*xx/(2*sig2*sig2))
   return,comp1+comp2
end
