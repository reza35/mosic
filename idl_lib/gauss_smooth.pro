;; $Id: //depot/Release/IDL_81/idl/idldir/lib/gauss_smooth.pro#1 $
;;
;; Copyright (c) 2010-2011, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;+
;; Gauss_Smooth
;;
;; Purpose:
;;   This function implements the gaussian smoothing function
;;
;-

;;---------------------------------------------------------------------------
;; Create_gaussian
;;
;; Purpose:
;;   Create a guassian kernal for smoothing
;;
;; Parameters:
;;   SIGMA - sigma value
;;
;;   WIDTH - desired width of the kernel
;;
;;   DIMS - number of dimensions of the kernel, currently only 1 or 2
;;
;; Keywords:
;;   NONE
;;
FUNCTION create_gaussian, sigmaIn, widthIn, dims

  ;; if not specified create a 5x5 1-sigma kernel
  sigma = (n_elements(sigmaIn) eq 0) ? 1.0d : double(sigmaIn)
  width = (n_elements(widthIn) eq 0) ? 5 : widthIn
  if (n_elements(dims) eq 0) then dims = 2
  
  ;; Ensure width is an odd integer
  width = fix(width) > 1
  width += width/2*2 eq width
  if ((n_elements(width) eq 1) && (dims eq 2)) then width = [width,width]
  if ((n_elements(sigma) eq 1) && (dims eq 2)) then sigma = [sigma,sigma]

  ;; create X and Y indices
  x = (dindgen(width[0])-width[0]/2)
  if (dims eq 2) then begin
    x #= replicate(1, width[1])
    y = transpose((dindgen(width[1])-width[1]/2) # $
      replicate(1, width[0]))
  endif

  ;; create kernel
  if (dims eq 2) then begin
    kernel = exp(-((x^2)/(2*sigma[0]^2)+(y^2)/(2*sigma[1]^2)))
  endif else begin
    kernel = exp(-(x^2)/(2*sigma^2))
  endelse

  return, kernel

END


;;---------------------------------------------------------------------------
;; Gauss_Smooth
;;
;; Purpose:
;;   Main routine
;;
FUNCTION gauss_smooth, data, sigmaIn, $
                       width=widthIn, kernel=gaussKernel, $
                       _EXTRA=_extra

on_error,2

  ;; verify we do have a 1D or 2D input
  dataDims = size(data, /n_dimensions)
  if ((dataDims lt 1) || (dataDims gt 2)) then begin
    message, 'Incorrect input data dimensions'
    return, 0
  endif
  ;; collapse dimensions of 1
  dims = size(data, /dimensions)
  dims1cnt = 0
  if (dataDims eq 2) then begin
    index = where(dims eq 1, dims1cnt)
    if (dims1cnt eq 1) then begin
      dataDims = 1
      ; Ensure no 1's are in any of the data dimensions
      data = reform(data, /OVERWRITE)
    endif
  endif
  
  ;; weed out inappropriate input types
  tname = size(data, /tname)
  if ((tname eq 'STRING') || (tname eq 'STRUCT') || (tname eq 'POINTER') || $
      (tname eq 'OBJREF')) then begin
    message, 'Incorrect input data type'
    return, 0
  endif

  ;; verify sigma
  sigmaDims = n_elements(sigmaIn)
  if (sigmaDims gt 2) then begin
    message, 'Incorrect sigma specification'
    return, 0
  endif
  if (sigmaDims eq 0) then $
    sigmaIn = 1.0d
  sigma = sigmaIn
  
  ;; weed out inappropriate sigma types
  tname = size(sigma, /tname)
  if ((tname eq 'STRING') || (tname eq 'STRUCT') || (tname eq 'POINTER') || $
      (tname eq 'OBJREF')) then begin
    message, 'Incorrect sigma type'
    return, 0
  endif

  if (sigmaDims gt dataDims) then begin
    message, 'Incompatible dimensions for Data and Sigma'
    return, 0
  endif
  
  ;; verify width
  widthDims = n_elements(widthIn)
  if (widthDims gt 2) then begin
    message, 'Incorrect width specification'
    return, 0
  endif
  ;; fill in width if not provided
  if (widthDims eq 0) then begin
    width = CEIL(sigma*3)
    ;; ensure width is odd
    width += width/2*2 eq width
    width *= 2
    width++
  endif else begin
    width = widthIn
  endelse
  
  ;; weed out inappropriate width types
  tname = size(widthIn, /tname)
  if ((tname eq 'STRING') || (tname eq 'STRUCT') || (tname eq 'POINTER') || $
      (tname eq 'OBJREF')) then begin
    message, 'Incorrect width type'
    return, 0
  endif

  if (widthDims gt dataDims) then begin
    message, 'Incompatible dimensions for Data and Width'
    return, 0
  endif
  
  if (max(size(data, /dimensions) lt width) ne 0) then begin
    message, 'Kernel size cannot be larger than the data size'
    return, 0
  endif

  ;; Ensure inputs fall within realistic values
  ((sigma >= 0.1)) <= 99.9

  dataType = size(data, /TYPE)
  gaussKernel = create_gaussian(sigma, width, dataDims)
  ; Kernel must be same data type as input
  gaussKernel = fix(temporary(gaussKernel), TYPE=dataType)

  smoothData = convol(data, gaussKernel, total(gaussKernel), $
                      _EXTRA=_extra)

  ; Put original data dimensions back                      
  if (dims1cnt eq 1) then $
    data = reform(data, dims, /OVERWRITE)

  return, smoothData

END
