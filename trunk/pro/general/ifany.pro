function ifany, arrcon

;Checks where any element from an input conditional is satisfied.  Returns 1 if any element satisfies the condition, 0 otherwise. Useful as a conditional on a loop.

;INPUTS
; arrcon - the conditonal array, of size N elements

;EXAMPLE
; x=[1,5,3,2,6]
; print, ifany(x mod 5 eq 0)
; returns, 1

; x=[1,5,3,2,6]
; if ifany(x ge 3) then begin loop

x=where(arrcon eq 1)
if x[0] ne -1 and n_elements(x) ge 1 then return, 1 else return, 0

end


function ifall, arrcon

;Checks where all elements from an input conditional are satisfied.  Returns 1 if all satisfy the condition, 0 otherwise. 

;INPUTS
;arrcon - the conditional array, of size N elements

return,  total(arrcon) eq n_elements(arrcon)

end