; Function allows new tags to be inserted into structure arrays.  Output is a new structure array
; with all variables copied and new empty fields for the new tags.
;
; KBundy 4/4/07
;
; SYNTAX:
;   new_str = adj_struct(old_struct, new_tag_names, new_default_values, [adjacent_tag_names],
;                        /BEFORE)
;     OR
;
;   new_str = adj_struct(old_struct, old_tag_names, new_tag_name, /RENAME)
;
;     OR
;
;   new_str = adj_struct(old_struct, delete_tag_names, /DELETE)
;
; INPUT
;   str:                The structure you want to insert tags into.
;   new_tag_names:      String array of new tag names.
;   new_default_values: Specifies out to fill the new empty tag fields.
;       If only inserting one new tag, this argument can be a scalar of any type.
;       If inserting more than one new tag, set this to a structure array with the same number of
;   elements as the number of new tags.
;
; OPTIONAL INPUT
;   adjacent_tag_names: If specified, new tags will be inserted before or after (the default) the
;                       original tags specified here.  Specify one for
;                       each new tag added.  Be sure to use upper case.
;
; KEYWORDS
;   BEFORE:  The default is to place new tags after those specified in adjacent_tag_names.  Setting
;            this keyword puts them before.
;   
; HISTORY
;   KBundy 4/10/07 - Added the RENAME keyword.  Can now rename tags in the structure.
;                    Added the DELETE keyword.  Can now delete tags.

function adj_struct, str, new_tag_names, new_default_values, adjacent_tag_names, $
                     BEFORE = BEFORE, RENAME = RENAME, DELETE = DELETE
  
; Warn if new_default_values is in pointer instead of structure format
; (old version was pointers)
if size(new_default_values, /type) NE 8 AND not(keyword_set(RENAME)) AND not(keyword_set(DELETE)) AND $
   n_elements(new_tag_names) NE 1 then begin
   print, '**FAILURE**  Program Updated: Please give new_default_values as a structure instead of a pointer'
   return, str
endif

n_str = n_elements(str)  ; Number of elements in structure array
str0 = str[0]   ; Work with the first element

; If we're just going to rename some of the tags...
if keyword_set(RENAME) then begin

    old_tag_names = new_tag_names      ; Rename the arguments (default arg names assume we're adding tags)
    new_tag_names = new_default_values ; Rename the arguments (default arg names assume we're adding tags)

    if n_elements(old_tag_names) NE n_elements(new_tag_names) then begin
        print, 'ERROR: Number of new and old tags should be the same.'
        return, -1
    endif

    n_tags_new = n_elements(new_tag_names)  ; Number of new tags being renamed
    n_tags_orig = n_tags(str0)        ; Original number of tags
    tag_names_orig = tag_names(str0)  ; Original tag names
    tag_names_output = tag_names_orig    ; Will contain updated tag name list
    
    w_rename = intarr(n_tags_new)
    for i = 0, n_tags_new-1 do w_rename[i] = where(old_tag_names[i] EQ tag_names_orig)
    w_badmatch = where(w_rename LT 0, n_badmatch)
    if n_badmatch GT 0 then begin
        print, 'ERROR: Some of the selected tag names are not contained in the input structure.'
        return, -1
    endif
    tag_names_output[w_rename] = new_tag_names  ; Replace with new tag names

    ; Determine the default values
    default_values = ptrarr(n_tags_orig, /allocate)
    for i = 0, n_tags_orig-1 do *default_values[i] = str0.(i)

    ; Create new 1st-element structure
    str0_output = create_struct(tag_names_output[0], *default_values[0])
    for i = 1, n_tags_orig-1 do str0_output = create_struct(str0_output, tag_names_output[i], *default_values[i])

    ; Make new structure array and fill with original values
    str_output = replicate(str0_output, n_str)
    for i = 0, n_tags_orig-1 do str_output.(i) = str.(i)

    return, str_output
endif

; If we're going to delete some of the tags
if keyword_set(DELETE) then begin
    tags_to_delete = new_tag_names  ; Rename the arguments (default arg names assume we're adding tags)
    n_delete = n_elements(tags_to_delete)

    n_tags_orig = n_tags(str0)        ; Original number of tags
    tag_names_orig = tag_names(str0)  ; Original tag names
    w_tag_orig = indgen(n_tags_orig)  ; Original index of the tags
    tag_names_output = tag_names_orig

    w_delete = intarr(n_delete)
    for i = 0, n_delete-1 do w_delete[i] = where(tags_to_delete[i] EQ tag_names_orig)
    w_badmatch = where(w_delete LT 0, n_badmatch)
    if n_badmatch GT 0 then begin
        print, 'ERROR: Some of the selected tag names are not contained in the input structure.'
        return, -1
    endif

    remove, w_delete, tag_names_output, w_tag_orig
    n_tags_output = n_elements(tag_names_output)

    ; Determine the default values
    default_values = ptrarr(n_tags_output, /allocate)
    for i = 0, n_tags_output-1 do *default_values[i] = str0.(w_tag_orig[i])

    ; Create new 1st-element structure
    str0_output = create_struct(tag_names_output[0], *default_values[0])
    for i = 1, n_tags_output-1 do str0_output = create_struct(str0_output, tag_names_output[i], *default_values[i])

    ; Make new structure array and fill with original values
    str_output = replicate(str0_output, n_str)
    for i = 0, n_tags_output-1 do str_output.(i) = str.(w_tag_orig[i])

    return, str_output
endif

n_tags_new = n_elements(new_tag_names)  ; Number of new tags being added
n_tags_orig = n_tags(str0)        ; Original number of tags
tag_names_orig = tag_names(str0)  ; Original tag names

wpos_mark = intarr(n_tags_new)  ; Tells you which of original tags the new ones should be next two
added_new = intarr(n_tags_new + n_tags_orig)  ; =1 where slot is to be filled, =0 everywhere else
tag_names_output = strarr(n_tags_new + n_tags_orig)  ; Outupt tag list with new ones added

; Determine where in the tag order the new tags should fall
if n_elements(adjacent_tag_names) GT 0 then begin
    for i = 0, n_tags_new-1 do wpos_mark[i] = where(tag_names_orig EQ adjacent_tag_names[i])

    ; This statement adjusts wpos_mark upward by the number of tags added
    ; before any given tag.  Only important if n_tags_new GT 1
    wpos_mark = wpos_mark + sort(wpos_mark)

;   If 'adjacent_tag_names' not specified, then place at the end of the list
endif else wpos_mark = indgen(n_tags_new) + (n_tags_orig-1)  

if keyword_set(BEFORE) then adjacent_position = 0 else adjacent_position = 1

added_new[wpos_mark+adjacent_position] = 1
w_keep_orig = where(added_new EQ 0)
w_new_tags = where(added_new EQ 1)

tag_names_output[w_keep_orig] = tag_names_orig
tag_names_output[w_new_tags] = new_tag_names

; Setup default values array
default_values = ptrarr(n_tags_new + n_tags_orig, /allocate)
for i = 0, n_tags_orig-1 do *default_values[w_keep_orig[i]] = str0.(i)

;if n_tags_new GT 1 then for i = 0, n_tags_new-1 do *default_values[w_new_tags[i]] = *new_default_values[i]  ; Assume new_default_values is a pointer here
if n_tags_new GT 1 then for i = 0, n_tags_new-1 do *default_values[w_new_tags[i]] = new_default_values.(i)  ; Assume new_default_values is a structure here
if n_tags_new EQ 1 then *default_values[w_new_tags[0]] = new_default_values                                  ; Assume new_default_values is not a pointer

; Create new 1st-element structure
str0_output = create_struct(tag_names_output[0], *default_values[0])
for i = 1, n_tags_new+n_tags_orig-1 do str0_output = create_struct(str0_output, tag_names_output[i], *default_values[i])

; Make new structure array and fill with original values
str_output = replicate(str0_output, n_str)
for i = 0, n_tags_orig-1 do str_output.(w_keep_orig[i]) = str.(i)

return, str_output

end
