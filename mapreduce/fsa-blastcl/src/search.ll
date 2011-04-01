; ModuleID = 'search.cl'
target datalayout = "e-p:32:32:32-i1:8:8-i8:8:8-i16:16:16-i32:32:32-i64:64:64-f32:32:32-f64:64:64-f80:32:32-v16:16:16-v24:32:32-v32:32:32-v48:64:64-v64:64:64-v96:128:128-v128:128:128-v192:256:256-v256:256:256-v512:512:512-v1024:1024:1024-a0:0:64"
target triple = "amdil-pc-amdopencl"

%0 = type { i8*, i8*, i8*, i8*, i32 }
%struct.PSSMatrixFP = type <{ <4 x i32>, <4 x i32>, <4 x i32>, <4 x i32>, <2 x i32>*, i8*, i8*, i8*, i8*, i32, i32, i32 }>
%struct.coordinate = type <{ <4 x i32>, <4 x i32> }>
%struct.groupFP = type <{ <4 x i32>, <4 x i32> }>
%struct.parameters = type <{ i8, i8, i8, i8, i8, i8, i8, i8, i8, i8, i8, i8, i8, i8, i8, i8, <4 x i32>, <4 x i32>, <4 x i32>, <4 x i32>, <4 x i32>, <4 x i32>, <4 x i32> }>
%struct.sequenceDataFP = type <{ <4 x i32>, <4 x i32>, <4 x i32>, <4 x i32>, <4 x i32> }>
%struct.ungappedExtension = type <{ %struct.coordinate, %struct.coordinate, %struct.coordinate, <4 x i32>, <4 x i32>, <4 x i32>, i8, i8, i8, i8, %struct.ungappedExtension*, i32, i32 }>

@sgv = internal constant [1 x i8] zeroinitializer ; <[1 x i8]*> [#uses=1]
@fgv = internal constant [1 x i8] zeroinitializer ; <[1 x i8]*> [#uses=1]
@lvgv = internal constant [0 x i8*] zeroinitializer ; <[0 x i8*]*> [#uses=1]
@global_numAdditionalTriggerExtensions = global <4 x i32> zeroinitializer, align 16 ; <<4 x i32>*> [#uses=1]
@sgv1 = internal constant [1 x i8] zeroinitializer ; <[1 x i8]*> [#uses=1]
@fgv2 = internal constant [1 x i8] zeroinitializer ; <[1 x i8]*> [#uses=1]
@lvgv3 = internal constant [0 x i8*] zeroinitializer ; <[0 x i8*]*> [#uses=1]
@llvm.global.annotations = appending global [2 x %0] [%0 { i8* bitcast (void (%struct.PSSMatrixFP addrspace(1)*, <2 x i32> addrspace(1)*, %struct.sequenceDataFP addrspace(1)*, i8 addrspace(1)*, %struct.parameters addrspace(1)*, %struct.coordinate addrspace(1)*, i8 addrspace(1)*, <4 x i32> addrspace(1)*, <4 x i32> addrspace(1)*, <4 x i32> addrspace(1)*, <4 x i32> addrspace(1)*, <4 x i32> addrspace(1)*, %struct.ungappedExtension addrspace(1)*, <4 x i32>)* @__OpenCL_search_protein1hitKernel_kernel to i8*), i8* getelementptr inbounds ([1 x i8]* @sgv, i32 0, i32 0), i8* getelementptr inbounds ([1 x i8]* @fgv, i32 0, i32 0), i8* bitcast ([0 x i8*]* @lvgv to i8*), i32 0 }, %0 { i8* bitcast (void (%struct.PSSMatrixFP addrspace(1)*, <2 x i32> addrspace(1)*, %struct.sequenceDataFP addrspace(1)*, i8 addrspace(1)*, %struct.parameters addrspace(1)*, %struct.coordinate addrspace(1)*, i8 addrspace(1)*, <4 x i32> addrspace(1)*, <4 x i32> addrspace(1)*, <4 x i32> addrspace(1)*, <4 x i32> addrspace(1)*, <4 x i32> addrspace(1)*, %struct.ungappedExtension addrspace(1)*, <4 x i32>)* @__OpenCL_search_protein2hitKernel_kernel to i8*), i8* getelementptr inbounds ([1 x i8]* @sgv1, i32 0, i32 0), i8* getelementptr inbounds ([1 x i8]* @fgv2, i32 0, i32 0), i8* bitcast ([0 x i8*]* @lvgv3 to i8*), i32 0 }], section "llvm.metadata" ; <[2 x %0]*> [#uses=0]

define void @__OpenCL_search_protein1hitKernel_kernel(%struct.PSSMatrixFP addrspace(1)* %PSSMatrixFP, <2 x i32> addrspace(1)* %matrixBody, %struct.sequenceDataFP addrspace(1)* %sequenceDataFP, i8 addrspace(1)* %sequence, %struct.parameters addrspace(1)* %parametersFP, %struct.coordinate addrspace(1)* %wordLookupDFA_groupsFP, i8 addrspace(1)* %wordLookupDFAFP, <4 x i32> addrspace(1)* %blast_numUngappedExtensions, <4 x i32> addrspace(1)* %blast_numTriggerExtensions, <4 x i32> addrspace(1)* %blast_numHits, <4 x i32> addrspace(1)* %hitMatrix_furthestp, <4 x i32> addrspace(1)* %hitMatrix_offsetp, %struct.ungappedExtension addrspace(1)* %ungappedExtension_extensionsp, <4 x i32> %nTotalSequenceNum) nounwind {
entry:
  %PSSMatrixFP.addr = alloca %struct.PSSMatrixFP addrspace(1)*, align 4 ; <%struct.PSSMatrixFP addrspace(1)**> [#uses=4]
  %matrixBody.addr = alloca <2 x i32> addrspace(1)*, align 4 ; <<2 x i32> addrspace(1)**> [#uses=2]
  %sequenceDataFP.addr = alloca %struct.sequenceDataFP addrspace(1)*, align 4 ; <%struct.sequenceDataFP addrspace(1)**> [#uses=4]
  %sequence.addr = alloca i8 addrspace(1)*, align 4 ; <i8 addrspace(1)**> [#uses=4]
  %parametersFP.addr = alloca %struct.parameters addrspace(1)*, align 4 ; <%struct.parameters addrspace(1)**> [#uses=10]
  %wordLookupDFA_groupsFP.addr = alloca %struct.coordinate addrspace(1)*, align 4 ; <%struct.coordinate addrspace(1)**> [#uses=6]
  %wordLookupDFAFP.addr = alloca i8 addrspace(1)*, align 4 ; <i8 addrspace(1)**> [#uses=3]
  %blast_numUngappedExtensions.addr = alloca <4 x i32> addrspace(1)*, align 4 ; <<4 x i32> addrspace(1)**> [#uses=2]
  %blast_numTriggerExtensions.addr = alloca <4 x i32> addrspace(1)*, align 4 ; <<4 x i32> addrspace(1)**> [#uses=2]
  %blast_numHits.addr = alloca <4 x i32> addrspace(1)*, align 4 ; <<4 x i32> addrspace(1)**> [#uses=2]
  %hitMatrix_furthestp.addr = alloca <4 x i32> addrspace(1)*, align 4 ; <<4 x i32> addrspace(1)**> [#uses=2]
  %hitMatrix_offsetp.addr = alloca <4 x i32> addrspace(1)*, align 4 ; <<4 x i32> addrspace(1)**> [#uses=2]
  %ungappedExtension_extensionsp.addr = alloca %struct.ungappedExtension addrspace(1)*, align 4 ; <%struct.ungappedExtension addrspace(1)**> [#uses=3]
  %nTotalSequenceNum.addr = alloca <4 x i32>, align 16 ; <<4 x i32>*> [#uses=2]
  %subject = alloca i8*, align 4                  ; <i8**> [#uses=4]
  %sequenceEnd = alloca i8*, align 4              ; <i8**> [#uses=3]
  %address = alloca i8*, align 4                  ; <i8**> [#uses=15]
  %subjectOffset = alloca <4 x i32>, align 16     ; <<4 x i32>*> [#uses=2]
  %count = alloca <4 x i32>, align 16             ; <<4 x i32>*> [#uses=4]
  %currentWord = alloca i8, align 1               ; <i8*> [#uses=4]
  %currentBlock = alloca i8*, align 4             ; <i8**> [#uses=4]
  %currentGroupFP = alloca %struct.coordinate*, align 4 ; <%struct.coordinate**> [#uses=10]
  %wordLookupDFA_AddiPositions = alloca <2 x i32>*, align 4 ; <<2 x i32>**> [#uses=2]
  %queryOffsets = alloca <2 x i32>*, align 4      ; <<2 x i32>**> [#uses=9]
  %queryOffset = alloca <2 x i32>, align 8        ; <<2 x i32>*> [#uses=3]
  %ungappedExtension_current = alloca %struct.ungappedExtension*, align 4 ; <%struct.ungappedExtension**> [#uses=2]
  %diagonal = alloca <4 x i32>, align 16          ; <<4 x i32>*> [#uses=2]
  %lastHitFP = alloca <4 x i32>*, align 4         ; <<4 x i32>**> [#uses=3]
  %ungappedExtension_subjectEndReachedFP = alloca <4 x i32>, align 16 ; <<4 x i32>*> [#uses=2]
  %hitMatrix_Local = alloca <4 x i32>*, align 4   ; <<4 x i32>**> [#uses=2]
  %sequenceCount = alloca <4 x i32>, align 16     ; <<4 x i32>*> [#uses=8]
  %bid = alloca i32, align 4                      ; <i32*> [#uses=2]
  %tid = alloca i32, align 4                      ; <i32*> [#uses=8]
  %numOfTriggerExtensions = alloca <4 x i32>, align 16 ; <<4 x i32>*> [#uses=3]
  %tmp224 = alloca <2 x i32>, align 8             ; <<2 x i32>*> [#uses=2]
  %tmp278 = alloca %struct.PSSMatrixFP, align 16  ; <%struct.PSSMatrixFP*> [#uses=2]
  store %struct.PSSMatrixFP addrspace(1)* %PSSMatrixFP, %struct.PSSMatrixFP addrspace(1)** %PSSMatrixFP.addr
  store <2 x i32> addrspace(1)* %matrixBody, <2 x i32> addrspace(1)** %matrixBody.addr
  store %struct.sequenceDataFP addrspace(1)* %sequenceDataFP, %struct.sequenceDataFP addrspace(1)** %sequenceDataFP.addr
  store i8 addrspace(1)* %sequence, i8 addrspace(1)** %sequence.addr
  store %struct.parameters addrspace(1)* %parametersFP, %struct.parameters addrspace(1)** %parametersFP.addr
  store %struct.coordinate addrspace(1)* %wordLookupDFA_groupsFP, %struct.coordinate addrspace(1)** %wordLookupDFA_groupsFP.addr
  store i8 addrspace(1)* %wordLookupDFAFP, i8 addrspace(1)** %wordLookupDFAFP.addr
  store <4 x i32> addrspace(1)* %blast_numUngappedExtensions, <4 x i32> addrspace(1)** %blast_numUngappedExtensions.addr
  store <4 x i32> addrspace(1)* %blast_numTriggerExtensions, <4 x i32> addrspace(1)** %blast_numTriggerExtensions.addr
  store <4 x i32> addrspace(1)* %blast_numHits, <4 x i32> addrspace(1)** %blast_numHits.addr
  store <4 x i32> addrspace(1)* %hitMatrix_furthestp, <4 x i32> addrspace(1)** %hitMatrix_furthestp.addr
  store <4 x i32> addrspace(1)* %hitMatrix_offsetp, <4 x i32> addrspace(1)** %hitMatrix_offsetp.addr
  store %struct.ungappedExtension addrspace(1)* %ungappedExtension_extensionsp, %struct.ungappedExtension addrspace(1)** %ungappedExtension_extensionsp.addr
  store <4 x i32> %nTotalSequenceNum, <4 x i32>* %nTotalSequenceNum.addr
  %call = call i32 @get_local_id(i32 1) nounwind  ; <i32> [#uses=1]
  %call1 = call i32 @get_global_size(i32 2) nounwind ; <i32> [#uses=1]
  %tmp = mul i32 %call, %call1                    ; <i32> [#uses=1]
  %call2 = call i32 @get_local_id(i32 2) nounwind ; <i32> [#uses=1]
  %tmp3 = add i32 %tmp, %call2                    ; <i32> [#uses=1]
  store i32 %tmp3, i32* %bid
  %tmp4 = load i32* %bid                          ; <i32> [#uses=1]
  %call5 = call i32 @get_local_size(i32 1) nounwind ; <i32> [#uses=1]
  %tmp6 = mul i32 %tmp4, %call5                   ; <i32> [#uses=1]
  %call7 = call i32 @get_local_size(i32 2) nounwind ; <i32> [#uses=1]
  %tmp8 = mul i32 %tmp6, %call7                   ; <i32> [#uses=1]
  %call9 = call i32 @get_group_id(i32 1) nounwind ; <i32> [#uses=1]
  %call10 = call i32 @get_local_size(i32 2) nounwind ; <i32> [#uses=1]
  %tmp11 = mul i32 %call9, %call10                ; <i32> [#uses=1]
  %tmp12 = add i32 %tmp8, %tmp11                  ; <i32> [#uses=1]
  %call13 = call i32 @get_group_id(i32 2) nounwind ; <i32> [#uses=1]
  %tmp14 = add i32 %tmp12, %call13                ; <i32> [#uses=1]
  store i32 %tmp14, i32* %tid
  store <4 x i32> zeroinitializer, <4 x i32>* %numOfTriggerExtensions
  %tmp15 = load <4 x i32> addrspace(1)** %hitMatrix_furthestp.addr ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp16 = load <4 x i32> addrspace(1)** %hitMatrix_offsetp.addr ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp17 = load i32* %tid                         ; <i32> [#uses=1]
  %arrayidx = getelementptr <4 x i32> addrspace(1)* %tmp16, i32 %tmp17 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp18 = load <4 x i32> addrspace(1)* %arrayidx ; <<4 x i32>> [#uses=1]
  %tmp19 = extractelement <4 x i32> %tmp18, i32 0 ; <i32> [#uses=1]
  %tmp20 = getelementptr <4 x i32> addrspace(1)* %tmp15, i32 %tmp19 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp21 = load %struct.PSSMatrixFP addrspace(1)** %PSSMatrixFP.addr ; <%struct.PSSMatrixFP addrspace(1)*> [#uses=1]
  %structele = getelementptr inbounds %struct.PSSMatrixFP addrspace(1)* %tmp21, i32 0, i32 0 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp22 = load <4 x i32> addrspace(1)* %structele ; <<4 x i32>> [#uses=1]
  %tmp23 = extractelement <4 x i32> %tmp22, i32 0 ; <i32> [#uses=1]
  %tmp24 = getelementptr <4 x i32> addrspace(1)* %tmp20, i32 %tmp23 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %conv = ptrtoint <4 x i32> addrspace(1)* %tmp24 to i32 ; <i32> [#uses=1]
  %conv25 = inttoptr i32 %conv to <4 x i32>*      ; <<4 x i32>*> [#uses=1]
  store <4 x i32>* %conv25, <4 x i32>** %hitMatrix_Local
  %tmp26 = load %struct.ungappedExtension addrspace(1)** %ungappedExtension_extensionsp.addr ; <%struct.ungappedExtension addrspace(1)*> [#uses=1]
  %structele27 = getelementptr inbounds %struct.ungappedExtension addrspace(1)* %tmp26, i32 0, i32 0 ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %structele28 = getelementptr inbounds %struct.coordinate addrspace(1)* %structele27, i32 0, i32 1 ; <<4 x i32> addrspace(1)*> [#uses=1]
  store <4 x i32> zeroinitializer, <4 x i32> addrspace(1)* %structele28
  %tmp29 = load %struct.ungappedExtension addrspace(1)** %ungappedExtension_extensionsp.addr ; <%struct.ungappedExtension addrspace(1)*> [#uses=1]
  %tmp30 = load i32* %tid                         ; <i32> [#uses=1]
  %tmp31 = mul i32 %tmp30, 150                    ; <i32> [#uses=1]
  %tmp32 = getelementptr %struct.ungappedExtension addrspace(1)* %tmp29, i32 %tmp31 ; <%struct.ungappedExtension addrspace(1)*> [#uses=1]
  %conv33 = ptrtoint %struct.ungappedExtension addrspace(1)* %tmp32 to i32 ; <i32> [#uses=1]
  %conv34 = inttoptr i32 %conv33 to %struct.ungappedExtension* ; <%struct.ungappedExtension*> [#uses=1]
  store %struct.ungappedExtension* %conv34, %struct.ungappedExtension** %ungappedExtension_current
  %tmp35 = load i8 addrspace(1)** %wordLookupDFAFP.addr ; <i8 addrspace(1)*> [#uses=1]
  %tmp36 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele37 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp36, i32 0, i32 17 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp38 = load <4 x i32> addrspace(1)* %structele37 ; <<4 x i32>> [#uses=1]
  %tmp39 = extractelement <4 x i32> %tmp38, i32 0 ; <i32> [#uses=1]
  %tmp40 = getelementptr i8 addrspace(1)* %tmp35, i32 %tmp39 ; <i8 addrspace(1)*> [#uses=1]
  %conv41 = ptrtoint i8 addrspace(1)* %tmp40 to i32 ; <i32> [#uses=1]
  %conv42 = inttoptr i32 %conv41 to <2 x i32>*    ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %conv42, <2 x i32>** %wordLookupDFA_AddiPositions
  %tmp43 = load %struct.PSSMatrixFP addrspace(1)** %PSSMatrixFP.addr ; <%struct.PSSMatrixFP addrspace(1)*> [#uses=1]
  %structele44 = getelementptr inbounds %struct.PSSMatrixFP addrspace(1)* %tmp43, i32 0, i32 4 ; <<2 x i32>* addrspace(1)*> [#uses=1]
  %tmp45 = load <2 x i32> addrspace(1)** %matrixBody.addr ; <<2 x i32> addrspace(1)*> [#uses=1]
  %tmp46 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele47 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp46, i32 0, i32 1 ; <i8 addrspace(1)*> [#uses=1]
  %tmp48 = load i8 addrspace(1)* %structele47     ; <i8> [#uses=1]
  %conv49 = zext i8 %tmp48 to i32                 ; <i32> [#uses=1]
  %tmp50 = getelementptr <2 x i32> addrspace(1)* %tmp45, i32 %conv49 ; <<2 x i32> addrspace(1)*> [#uses=1]
  %conv51 = ptrtoint <2 x i32> addrspace(1)* %tmp50 to i32 ; <i32> [#uses=1]
  %conv52 = inttoptr i32 %conv51 to <2 x i32>*    ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %conv52, <2 x i32>* addrspace(1)* %structele44
  %tmp53 = load i32* %tid                         ; <i32> [#uses=4]
  %conv54 = insertelement <4 x i32> undef, i32 %tmp53, i32 0 ; <<4 x i32>> [#uses=1]
  %conv55 = insertelement <4 x i32> %conv54, i32 %tmp53, i32 1 ; <<4 x i32>> [#uses=1]
  %conv56 = insertelement <4 x i32> %conv55, i32 %tmp53, i32 2 ; <<4 x i32>> [#uses=1]
  %conv57 = insertelement <4 x i32> %conv56, i32 %tmp53, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv57, <4 x i32>* %sequenceCount
  br label %while.cond

return:                                           ; preds = %while.exit
  ret void

while.cond:                                       ; preds = %if.end, %entry
  %tmp58 = load <4 x i32>* %sequenceCount         ; <<4 x i32>> [#uses=1]
  %tmp59 = extractelement <4 x i32> %tmp58, i32 0 ; <i32> [#uses=1]
  %tmp60 = load <4 x i32>* %nTotalSequenceNum.addr ; <<4 x i32>> [#uses=1]
  %tmp61 = extractelement <4 x i32> %tmp60, i32 0 ; <i32> [#uses=1]
  %cmp = icmp ult i32 %tmp59, %tmp61              ; <i1> [#uses=1]
  br i1 %cmp, label %while.body, label %while.exit

while.exit:                                       ; preds = %while.cond
  %tmp321 = load <4 x i32> addrspace(1)** %blast_numTriggerExtensions.addr ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp322 = load i32* %tid                        ; <i32> [#uses=1]
  %arrayidx323 = getelementptr <4 x i32> addrspace(1)* %tmp321, i32 %tmp322 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp324 = load <4 x i32>* %numOfTriggerExtensions ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp324, <4 x i32> addrspace(1)* %arrayidx323
  br label %return

while.body:                                       ; preds = %while.cond
  %tmp62 = load i8 addrspace(1)** %sequence.addr  ; <i8 addrspace(1)*> [#uses=1]
  %tmp63 = load %struct.sequenceDataFP addrspace(1)** %sequenceDataFP.addr ; <%struct.sequenceDataFP addrspace(1)*> [#uses=1]
  %tmp64 = load <4 x i32>* %sequenceCount         ; <<4 x i32>> [#uses=1]
  %tmp65 = extractelement <4 x i32> %tmp64, i32 0 ; <i32> [#uses=1]
  %arrayidx66 = getelementptr %struct.sequenceDataFP addrspace(1)* %tmp63, i32 %tmp65 ; <%struct.sequenceDataFP addrspace(1)*> [#uses=1]
  %structele67 = getelementptr inbounds %struct.sequenceDataFP addrspace(1)* %arrayidx66, i32 0, i32 4 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp68 = load <4 x i32> addrspace(1)* %structele67 ; <<4 x i32>> [#uses=1]
  %tmp69 = extractelement <4 x i32> %tmp68, i32 0 ; <i32> [#uses=1]
  %tmp70 = getelementptr i8 addrspace(1)* %tmp62, i32 %tmp69 ; <i8 addrspace(1)*> [#uses=1]
  %conv71 = ptrtoint i8 addrspace(1)* %tmp70 to i32 ; <i32> [#uses=1]
  %conv72 = inttoptr i32 %conv71 to i8*           ; <i8*> [#uses=2]
  store i8* %conv72, i8** %address
  store i8* %conv72, i8** %subject
  %tmp73 = load %struct.sequenceDataFP addrspace(1)** %sequenceDataFP.addr ; <%struct.sequenceDataFP addrspace(1)*> [#uses=1]
  %tmp74 = load <4 x i32>* %sequenceCount         ; <<4 x i32>> [#uses=1]
  %tmp75 = extractelement <4 x i32> %tmp74, i32 0 ; <i32> [#uses=1]
  %arrayidx76 = getelementptr %struct.sequenceDataFP addrspace(1)* %tmp73, i32 %tmp75 ; <%struct.sequenceDataFP addrspace(1)*> [#uses=1]
  %structele77 = getelementptr inbounds %struct.sequenceDataFP addrspace(1)* %arrayidx76, i32 0, i32 0 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp78 = load <4 x i32> addrspace(1)* %structele77 ; <<4 x i32>> [#uses=1]
  %tmp79 = extractelement <4 x i32> %tmp78, i32 0 ; <i32> [#uses=1]
  %tmp80 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele81 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp80, i32 0, i32 0 ; <i8 addrspace(1)*> [#uses=1]
  %tmp82 = load i8 addrspace(1)* %structele81     ; <i8> [#uses=1]
  %conv83 = sext i8 %tmp82 to i32                 ; <i32> [#uses=1]
  %cmp84 = icmp uge i32 %tmp79, %conv83           ; <i1> [#uses=1]
  br i1 %cmp84, label %if.then, label %if.end

if.end:                                           ; preds = %__T148432052, %while.body
  %tmp312 = load <4 x i32>* %sequenceCount        ; <<4 x i32>> [#uses=1]
  %call313 = call i32 @get_global_size(i32 1) nounwind ; <i32> [#uses=1]
  %call314 = call i32 @get_local_size(i32 1) nounwind ; <i32> [#uses=1]
  %tmp315 = mul i32 %call313, %call314            ; <i32> [#uses=4]
  %conv316 = insertelement <4 x i32> undef, i32 %tmp315, i32 0 ; <<4 x i32>> [#uses=1]
  %conv317 = insertelement <4 x i32> %conv316, i32 %tmp315, i32 1 ; <<4 x i32>> [#uses=1]
  %conv318 = insertelement <4 x i32> %conv317, i32 %tmp315, i32 2 ; <<4 x i32>> [#uses=1]
  %conv319 = insertelement <4 x i32> %conv318, i32 %tmp315, i32 3 ; <<4 x i32>> [#uses=1]
  %tmp320 = add <4 x i32> %tmp312, %conv319       ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp320, <4 x i32>* %sequenceCount
  br label %while.cond

if.then:                                          ; preds = %while.body
  %tmp85 = load %struct.coordinate addrspace(1)** %wordLookupDFA_groupsFP.addr ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %conv86 = ptrtoint %struct.coordinate addrspace(1)* %tmp85 to i32 ; <i32> [#uses=1]
  %conv87 = inttoptr i32 %conv86 to %struct.coordinate* ; <%struct.coordinate*> [#uses=1]
  store %struct.coordinate* %conv87, %struct.coordinate** %currentGroupFP
  store <4 x i32> <i32 1, i32 1, i32 1, i32 1>, <4 x i32>* %count
  br label %while.cond88

while.cond88:                                     ; preds = %if.end106, %if.then
  %tmp89 = load <4 x i32>* %count                 ; <<4 x i32>> [#uses=1]
  %tmp90 = extractelement <4 x i32> %tmp89, i32 0 ; <i32> [#uses=1]
  %tmp91 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele92 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp91, i32 0, i32 0 ; <i8 addrspace(1)*> [#uses=1]
  %tmp93 = load i8 addrspace(1)* %structele92     ; <i8> [#uses=1]
  %conv94 = sext i8 %tmp93 to i32                 ; <i32> [#uses=1]
  %cmp95 = icmp slt i32 %tmp90, %conv94           ; <i1> [#uses=1]
  br i1 %cmp95, label %while.body97, label %while.exit96

while.exit96:                                     ; preds = %while.cond88
  %tmp132 = load i8** %subject                    ; <i8*> [#uses=1]
  %tmp133 = load %struct.sequenceDataFP addrspace(1)** %sequenceDataFP.addr ; <%struct.sequenceDataFP addrspace(1)*> [#uses=1]
  %tmp134 = load <4 x i32>* %sequenceCount        ; <<4 x i32>> [#uses=1]
  %tmp135 = extractelement <4 x i32> %tmp134, i32 0 ; <i32> [#uses=1]
  %arrayidx136 = getelementptr %struct.sequenceDataFP addrspace(1)* %tmp133, i32 %tmp135 ; <%struct.sequenceDataFP addrspace(1)*> [#uses=1]
  %structele137 = getelementptr inbounds %struct.sequenceDataFP addrspace(1)* %arrayidx136, i32 0, i32 0 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp138 = load <4 x i32> addrspace(1)* %structele137 ; <<4 x i32>> [#uses=1]
  %tmp139 = extractelement <4 x i32> %tmp138, i32 0 ; <i32> [#uses=1]
  %tmp140 = getelementptr i8* %tmp132, i32 %tmp139 ; <i8*> [#uses=1]
  store i8* %tmp140, i8** %sequenceEnd
  br label %while.cond141

while.body97:                                     ; preds = %while.cond88
  %tmp98 = load i8** %address                     ; <i8*> [#uses=1]
  %tmp99 = load i8* %tmp98                        ; <i8> [#uses=1]
  %conv100 = zext i8 %tmp99 to i32                ; <i32> [#uses=1]
  %tmp101 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele102 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp101, i32 0, i32 16 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp103 = load <4 x i32> addrspace(1)* %structele102 ; <<4 x i32>> [#uses=1]
  %tmp104 = extractelement <4 x i32> %tmp103, i32 0 ; <i32> [#uses=1]
  %cmp105 = icmp slt i32 %conv100, %tmp104        ; <i1> [#uses=1]
  br i1 %cmp105, label %if.then107, label %if.else

if.end106:                                        ; preds = %if.else, %if.then107
  %tmp128 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp129 = getelementptr i8* %tmp128, i32 1      ; <i8*> [#uses=1]
  store i8* %tmp129, i8** %address
  %tmp130 = load <4 x i32>* %count                ; <<4 x i32>> [#uses=1]
  %tmp131 = add <4 x i32> %tmp130, <i32 1, i32 1, i32 1, i32 1> ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp131, <4 x i32>* %count
  br label %while.cond88

if.then107:                                       ; preds = %while.body97
  %tmp108 = load %struct.coordinate addrspace(1)** %wordLookupDFA_groupsFP.addr ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %tmp109 = load %struct.coordinate** %currentGroupFP ; <%struct.coordinate*> [#uses=1]
  %structele110 = getelementptr inbounds %struct.coordinate* %tmp109, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp111 = load <4 x i32>* %structele110         ; <<4 x i32>> [#uses=1]
  %tmp112 = extractelement <4 x i32> %tmp111, i32 0 ; <i32> [#uses=1]
  %tmp113 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp114 = load i8* %tmp113                      ; <i8> [#uses=1]
  %conv115 = zext i8 %tmp114 to i32               ; <i32> [#uses=1]
  %tmp116 = add i32 %tmp112, %conv115             ; <i32> [#uses=1]
  %tmp117 = getelementptr %struct.coordinate addrspace(1)* %tmp108, i32 %tmp116 ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %conv118 = ptrtoint %struct.coordinate addrspace(1)* %tmp117 to i32 ; <i32> [#uses=1]
  %conv119 = inttoptr i32 %conv118 to %struct.coordinate* ; <%struct.coordinate*> [#uses=1]
  store %struct.coordinate* %conv119, %struct.coordinate** %currentGroupFP
  br label %if.end106

if.else:                                          ; preds = %while.body97
  %tmp120 = load %struct.coordinate addrspace(1)** %wordLookupDFA_groupsFP.addr ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %tmp121 = load %struct.coordinate** %currentGroupFP ; <%struct.coordinate*> [#uses=1]
  %structele122 = getelementptr inbounds %struct.coordinate* %tmp121, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp123 = load <4 x i32>* %structele122         ; <<4 x i32>> [#uses=1]
  %tmp124 = extractelement <4 x i32> %tmp123, i32 0 ; <i32> [#uses=1]
  %tmp125 = getelementptr %struct.coordinate addrspace(1)* %tmp120, i32 %tmp124 ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %conv126 = ptrtoint %struct.coordinate addrspace(1)* %tmp125 to i32 ; <i32> [#uses=1]
  %conv127 = inttoptr i32 %conv126 to %struct.coordinate* ; <%struct.coordinate*> [#uses=1]
  store %struct.coordinate* %conv127, %struct.coordinate** %currentGroupFP
  br label %if.end106

while.cond141:                                    ; preds = %if.end200, %while.exit96
  %tmp142 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp143 = load i8** %sequenceEnd                ; <i8*> [#uses=1]
  %cmp144 = icmp ult i8* %tmp142, %tmp143         ; <i1> [#uses=1]
  br i1 %cmp144, label %while.body146, label %while.exit145

while.exit145:                                    ; preds = %while.cond141
  br label %__T148432052

while.body146:                                    ; preds = %while.cond141
  %tmp147 = load i8 addrspace(1)** %wordLookupDFAFP.addr ; <i8 addrspace(1)*> [#uses=1]
  %tmp148 = load %struct.coordinate** %currentGroupFP ; <%struct.coordinate*> [#uses=1]
  %structele149 = getelementptr inbounds %struct.coordinate* %tmp148, i32 0, i32 0 ; <<4 x i32>*> [#uses=1]
  %tmp150 = load <4 x i32>* %structele149         ; <<4 x i32>> [#uses=1]
  %tmp151 = extractelement <4 x i32> %tmp150, i32 0 ; <i32> [#uses=1]
  %tmp152 = getelementptr i8 addrspace(1)* %tmp147, i32 %tmp151 ; <i8 addrspace(1)*> [#uses=1]
  %conv153 = ptrtoint i8 addrspace(1)* %tmp152 to i32 ; <i32> [#uses=1]
  %conv154 = inttoptr i32 %conv153 to i8*         ; <i8*> [#uses=1]
  store i8* %conv154, i8** %currentBlock
  %tmp155 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp156 = load i8* %tmp155                      ; <i8> [#uses=1]
  %conv157 = zext i8 %tmp156 to i32               ; <i32> [#uses=1]
  %tmp158 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele159 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp158, i32 0, i32 16 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp160 = load <4 x i32> addrspace(1)* %structele159 ; <<4 x i32>> [#uses=1]
  %tmp161 = extractelement <4 x i32> %tmp160, i32 0 ; <i32> [#uses=1]
  %cmp162 = icmp slt i32 %conv157, %tmp161        ; <i1> [#uses=1]
  br i1 %cmp162, label %if.then164, label %if.else165

if.end163:                                        ; preds = %if.end187, %if.then164
  %tmp199 = load i8* %currentWord                 ; <i8> [#uses=1]
  %tobool = icmp ne i8 %tmp199, 0                 ; <i1> [#uses=1]
  br i1 %tobool, label %if.then201, label %if.end200

if.then164:                                       ; preds = %while.body146
  %tmp166 = load i8** %currentBlock               ; <i8*> [#uses=1]
  %tmp167 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp168 = load i8* %tmp167                      ; <i8> [#uses=1]
  %conv169 = zext i8 %tmp168 to i32               ; <i32> [#uses=1]
  %arrayidx170 = getelementptr i8* %tmp166, i32 %conv169 ; <i8*> [#uses=1]
  %tmp171 = load i8* %arrayidx170                 ; <i8> [#uses=1]
  store i8 %tmp171, i8* %currentWord
  %tmp172 = load %struct.coordinate addrspace(1)** %wordLookupDFA_groupsFP.addr ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %tmp173 = load %struct.coordinate** %currentGroupFP ; <%struct.coordinate*> [#uses=1]
  %structele174 = getelementptr inbounds %struct.coordinate* %tmp173, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp175 = load <4 x i32>* %structele174         ; <<4 x i32>> [#uses=1]
  %tmp176 = extractelement <4 x i32> %tmp175, i32 0 ; <i32> [#uses=1]
  %tmp177 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp178 = load i8* %tmp177                      ; <i8> [#uses=1]
  %conv179 = zext i8 %tmp178 to i32               ; <i32> [#uses=1]
  %tmp180 = add i32 %tmp176, %conv179             ; <i32> [#uses=1]
  %tmp181 = getelementptr %struct.coordinate addrspace(1)* %tmp172, i32 %tmp180 ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %conv182 = ptrtoint %struct.coordinate addrspace(1)* %tmp181 to i32 ; <i32> [#uses=1]
  %conv183 = inttoptr i32 %conv182 to %struct.coordinate* ; <%struct.coordinate*> [#uses=1]
  store %struct.coordinate* %conv183, %struct.coordinate** %currentGroupFP
  br label %if.end163

if.else165:                                       ; preds = %while.body146
  %tmp184 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp185 = load i8** %sequenceEnd                ; <i8*> [#uses=1]
  %cmp186 = icmp uge i8* %tmp184, %tmp185         ; <i1> [#uses=1]
  br i1 %cmp186, label %if.then188, label %if.end187

if.end187:                                        ; preds = %if.else165
  %tmp189 = load i8** %currentBlock               ; <i8*> [#uses=1]
  %tmp190 = load i8* %tmp189                      ; <i8> [#uses=1]
  store i8 %tmp190, i8* %currentWord
  %tmp191 = load %struct.coordinate addrspace(1)** %wordLookupDFA_groupsFP.addr ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %tmp192 = load %struct.coordinate** %currentGroupFP ; <%struct.coordinate*> [#uses=1]
  %structele193 = getelementptr inbounds %struct.coordinate* %tmp192, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp194 = load <4 x i32>* %structele193         ; <<4 x i32>> [#uses=1]
  %tmp195 = extractelement <4 x i32> %tmp194, i32 0 ; <i32> [#uses=1]
  %tmp196 = getelementptr %struct.coordinate addrspace(1)* %tmp191, i32 %tmp195 ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %conv197 = ptrtoint %struct.coordinate addrspace(1)* %tmp196 to i32 ; <i32> [#uses=1]
  %conv198 = inttoptr i32 %conv197 to %struct.coordinate* ; <%struct.coordinate*> [#uses=1]
  store %struct.coordinate* %conv198, %struct.coordinate** %currentGroupFP
  br label %if.end163

if.then188:                                       ; preds = %if.else165
  br label %__T148432052

__T148432052:                                     ; preds = %while.exit145, %if.then188
  br label %if.end

if.end200:                                        ; preds = %do.exit, %if.end163
  %tmp310 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp311 = getelementptr i8* %tmp310, i32 1      ; <i8*> [#uses=1]
  store i8* %tmp311, i8** %address
  br label %while.cond141

if.then201:                                       ; preds = %if.end163
  %tmp202 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp203 = load i8** %subject                    ; <i8*> [#uses=1]
  %sub.ptr.lhs.cast = ptrtoint i8* %tmp202 to i32 ; <i32> [#uses=1]
  %sub.ptr.rhs.cast = ptrtoint i8* %tmp203 to i32 ; <i32> [#uses=1]
  %sub.ptr.sub = sub i32 %sub.ptr.lhs.cast, %sub.ptr.rhs.cast ; <i32> [#uses=4]
  %conv204 = insertelement <4 x i32> undef, i32 %sub.ptr.sub, i32 0 ; <<4 x i32>> [#uses=1]
  %conv205 = insertelement <4 x i32> %conv204, i32 %sub.ptr.sub, i32 1 ; <<4 x i32>> [#uses=1]
  %conv206 = insertelement <4 x i32> %conv205, i32 %sub.ptr.sub, i32 2 ; <<4 x i32>> [#uses=1]
  %conv207 = insertelement <4 x i32> %conv206, i32 %sub.ptr.sub, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv207, <4 x i32>* %subjectOffset
  %tmp208 = load i8** %currentBlock               ; <i8*> [#uses=1]
  %conv209 = bitcast i8* %tmp208 to <2 x i32>*    ; <<2 x i32>*> [#uses=1]
  %tmp210 = load i8* %currentWord                 ; <i8> [#uses=1]
  %conv211 = zext i8 %tmp210 to i32               ; <i32> [#uses=1]
  %sub.ptr.neg = sub i32 0, %conv211              ; <i32> [#uses=1]
  %tmp212 = getelementptr <2 x i32>* %conv209, i32 %sub.ptr.neg ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp212, <2 x i32>** %queryOffsets
  %tmp213 = load <2 x i32>** %queryOffsets        ; <<2 x i32>*> [#uses=1]
  %tmp214 = load <2 x i32>* %tmp213               ; <<2 x i32>> [#uses=1]
  %tmp215 = extractelement <2 x i32> %tmp214, i32 0 ; <i32> [#uses=1]
  %tobool216 = icmp ne i32 %tmp215, 0             ; <i1> [#uses=1]
  %lnot = xor i1 %tobool216, true                 ; <i1> [#uses=1]
  br i1 %lnot, label %if.then218, label %if.end217

if.end217:                                        ; preds = %if.then218, %if.then201
  br label %do.body

if.then218:                                       ; preds = %if.then201
  %tmp219 = load <2 x i32>** %wordLookupDFA_AddiPositions ; <<2 x i32>*> [#uses=1]
  %tmp220 = load <2 x i32>** %queryOffsets        ; <<2 x i32>*> [#uses=1]
  %tmp221 = getelementptr <2 x i32>* %tmp220, i32 1 ; <<2 x i32>*> [#uses=1]
  %tmp222 = load <2 x i32>* %tmp221               ; <<2 x i32>> [#uses=1]
  %tmp223 = mul <2 x i32> %tmp222, <i32 32767, i32 32767> ; <<2 x i32>> [#uses=1]
  store <2 x i32> %tmp223, <2 x i32>* %tmp224
  %tmp225 = load <2 x i32>* %tmp224               ; <<2 x i32>> [#uses=1]
  %tmp226 = extractelement <2 x i32> %tmp225, i32 0 ; <i32> [#uses=1]
  %tmp227 = getelementptr <2 x i32>* %tmp219, i32 %tmp226 ; <<2 x i32>*> [#uses=1]
  %tmp228 = load <2 x i32>** %queryOffsets        ; <<2 x i32>*> [#uses=1]
  %tmp229 = getelementptr <2 x i32>* %tmp228, i32 2 ; <<2 x i32>*> [#uses=1]
  %tmp230 = load <2 x i32>* %tmp229               ; <<2 x i32>> [#uses=1]
  %tmp231 = extractelement <2 x i32> %tmp230, i32 0 ; <i32> [#uses=1]
  %tmp232 = getelementptr <2 x i32>* %tmp227, i32 %tmp231 ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp232, <2 x i32>** %queryOffsets
  br label %if.end217

do.exit:                                          ; preds = %do.cond
  br label %if.end200

do.body:                                          ; preds = %do.cond, %if.end217
  %tmp233 = load <2 x i32>** %queryOffsets        ; <<2 x i32>*> [#uses=1]
  %tmp234 = load <2 x i32>* %tmp233               ; <<2 x i32>> [#uses=1]
  store <2 x i32> %tmp234, <2 x i32>* %queryOffset
  %tmp235 = load <4 x i32>* %subjectOffset        ; <<4 x i32>> [#uses=1]
  %tmp236 = extractelement <4 x i32> %tmp235, i32 0 ; <i32> [#uses=1]
  %tmp237 = load <2 x i32>* %queryOffset          ; <<2 x i32>> [#uses=1]
  %tmp238 = extractelement <2 x i32> %tmp237, i32 0 ; <i32> [#uses=1]
  %tmp239 = sub i32 %tmp236, %tmp238              ; <i32> [#uses=4]
  %conv240 = insertelement <4 x i32> undef, i32 %tmp239, i32 0 ; <<4 x i32>> [#uses=1]
  %conv241 = insertelement <4 x i32> %conv240, i32 %tmp239, i32 1 ; <<4 x i32>> [#uses=1]
  %conv242 = insertelement <4 x i32> %conv241, i32 %tmp239, i32 2 ; <<4 x i32>> [#uses=1]
  %conv243 = insertelement <4 x i32> %conv242, i32 %tmp239, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv243, <4 x i32>* %diagonal
  %tmp244 = load <4 x i32>** %hitMatrix_Local     ; <<4 x i32>*> [#uses=1]
  %tmp245 = load <4 x i32>* %diagonal             ; <<4 x i32>> [#uses=1]
  %tmp246 = extractelement <4 x i32> %tmp245, i32 0 ; <i32> [#uses=1]
  %tmp247 = getelementptr <4 x i32>* %tmp244, i32 %tmp246 ; <<4 x i32>*> [#uses=1]
  store <4 x i32>* %tmp247, <4 x i32>** %lastHitFP
  %tmp248 = load <4 x i32>** %lastHitFP           ; <<4 x i32>*> [#uses=1]
  %tmp249 = load <4 x i32>* %tmp248               ; <<4 x i32>> [#uses=1]
  %tmp250 = extractelement <4 x i32> %tmp249, i32 0 ; <i32> [#uses=1]
  %tmp251 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp252 = load i8 addrspace(1)** %sequence.addr ; <i8 addrspace(1)*> [#uses=1]
  %conv253 = ptrtoint i8 addrspace(1)* %tmp252 to i32 ; <i32> [#uses=1]
  %conv254 = inttoptr i32 %conv253 to i8*         ; <i8*> [#uses=1]
  %sub.ptr.lhs.cast255 = ptrtoint i8* %tmp251 to i32 ; <i32> [#uses=1]
  %sub.ptr.rhs.cast256 = ptrtoint i8* %conv254 to i32 ; <i32> [#uses=1]
  %sub.ptr.sub257 = sub i32 %sub.ptr.lhs.cast255, %sub.ptr.rhs.cast256 ; <i32> [#uses=1]
  %cmp258 = icmp ult i32 %tmp250, %sub.ptr.sub257 ; <i1> [#uses=1]
  br i1 %cmp258, label %if.then260, label %if.end259

do.cond:                                          ; preds = %if.end259
  %tmp306 = load <2 x i32>** %queryOffsets        ; <<2 x i32>*> [#uses=1]
  %tmp307 = load <2 x i32>* %tmp306               ; <<2 x i32>> [#uses=1]
  %tmp308 = extractelement <2 x i32> %tmp307, i32 0 ; <i32> [#uses=1]
  %tobool309 = icmp ne i32 %tmp308, 0             ; <i1> [#uses=1]
  br i1 %tobool309, label %do.body, label %do.exit

if.end259:                                        ; preds = %if.then260, %do.body
  %tmp299 = load <2 x i32>** %queryOffsets        ; <<2 x i32>*> [#uses=1]
  %tmp300 = getelementptr <2 x i32>* %tmp299, i32 1 ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp300, <2 x i32>** %queryOffsets
  %tmp301 = load <4 x i32> addrspace(1)** %blast_numHits.addr ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp302 = load i32* %tid                        ; <i32> [#uses=1]
  %arrayidx303 = getelementptr <4 x i32> addrspace(1)* %tmp301, i32 %tmp302 ; <<4 x i32> addrspace(1)*> [#uses=2]
  %tmp304 = load <4 x i32> addrspace(1)* %arrayidx303 ; <<4 x i32>> [#uses=1]
  %tmp305 = add <4 x i32> %tmp304, <i32 1, i32 1, i32 1, i32 1> ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp305, <4 x i32> addrspace(1)* %arrayidx303
  br label %do.cond

if.then260:                                       ; preds = %do.body
  %tmp261 = load <4 x i32> addrspace(1)** %blast_numUngappedExtensions.addr ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp262 = load i32* %tid                        ; <i32> [#uses=1]
  %arrayidx263 = getelementptr <4 x i32> addrspace(1)* %tmp261, i32 %tmp262 ; <<4 x i32> addrspace(1)*> [#uses=2]
  %tmp264 = load <4 x i32> addrspace(1)* %arrayidx263 ; <<4 x i32>> [#uses=1]
  %tmp265 = add <4 x i32> %tmp264, <i32 1, i32 1, i32 1, i32 1> ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp265, <4 x i32> addrspace(1)* %arrayidx263
  %tmp266 = load i8 addrspace(1)** %sequence.addr ; <i8 addrspace(1)*> [#uses=1]
  %conv267 = ptrtoint i8 addrspace(1)* %tmp266 to i32 ; <i32> [#uses=1]
  %conv268 = inttoptr i32 %conv267 to i8*         ; <i8*> [#uses=1]
  %tmp269 = load <2 x i32>* %queryOffset          ; <<2 x i32>> [#uses=1]
  %tmp270 = extractelement <2 x i32> %tmp269, i32 0 ; <i32> [#uses=4]
  %conv271 = insertelement <4 x i32> undef, i32 %tmp270, i32 0 ; <<4 x i32>> [#uses=1]
  %conv272 = insertelement <4 x i32> %conv271, i32 %tmp270, i32 1 ; <<4 x i32>> [#uses=1]
  %conv273 = insertelement <4 x i32> %conv272, i32 %tmp270, i32 2 ; <<4 x i32>> [#uses=1]
  %conv274 = insertelement <4 x i32> %conv273, i32 %tmp270, i32 3 ; <<4 x i32>> [#uses=1]
  %tmp275 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp276 = load %struct.PSSMatrixFP addrspace(1)** %PSSMatrixFP.addr ; <%struct.PSSMatrixFP addrspace(1)*> [#uses=1]
  %tmp277 = load %struct.PSSMatrixFP addrspace(1)* %tmp276 ; <%struct.PSSMatrixFP> [#uses=1]
  store %struct.PSSMatrixFP %tmp277, %struct.PSSMatrixFP* %tmp278
  %tmp279 = load i8** %subject                    ; <i8*> [#uses=1]
  %tmp280 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele281 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp280, i32 0, i32 1 ; <i8 addrspace(1)*> [#uses=1]
  %tmp282 = load i8 addrspace(1)* %structele281   ; <i8> [#uses=1]
  %tmp283 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele284 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp283, i32 0, i32 18 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp285 = load <4 x i32> addrspace(1)* %structele284 ; <<4 x i32>> [#uses=1]
  %tmp286 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele287 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp286, i32 0, i32 19 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp288 = load <4 x i32> addrspace(1)* %structele287 ; <<4 x i32>> [#uses=1]
  %tmp289 = load %struct.ungappedExtension** %ungappedExtension_current ; <%struct.ungappedExtension*> [#uses=1]
  %tmp290 = load <4 x i32>* %sequenceCount        ; <<4 x i32>> [#uses=1]
  %tmp291 = load i32* %tid                        ; <i32> [#uses=4]
  %conv292 = insertelement <4 x i32> undef, i32 %tmp291, i32 0 ; <<4 x i32>> [#uses=1]
  %conv293 = insertelement <4 x i32> %conv292, i32 %tmp291, i32 1 ; <<4 x i32>> [#uses=1]
  %conv294 = insertelement <4 x i32> %conv293, i32 %tmp291, i32 2 ; <<4 x i32>> [#uses=1]
  %conv295 = insertelement <4 x i32> %conv294, i32 %tmp291, i32 3 ; <<4 x i32>> [#uses=1]
  %call296 = call %struct.ungappedExtension addrspace(1)* @ungappedExtension_oneHitExtendD(i8* %conv268, <4 x i32> %conv274, i8* %tmp275, %struct.PSSMatrixFP* byval %tmp278, i8* %tmp279, <4 x i32>* %ungappedExtension_subjectEndReachedFP, i8 zeroext %tmp282, <4 x i32> %tmp285, <4 x i32> %tmp288, %struct.ungappedExtension* %tmp289, <4 x i32>* %numOfTriggerExtensions, <4 x i32> %tmp290, <4 x i32> %conv295) nounwind ; <%struct.ungappedExtension addrspace(1)*> [#uses=0]
  %tmp297 = load <4 x i32>** %lastHitFP           ; <<4 x i32>*> [#uses=1]
  %tmp298 = load <4 x i32>* %ungappedExtension_subjectEndReachedFP ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp298, <4 x i32>* %tmp297
  br label %if.end259
}

declare i32 @get_local_id(i32) nounwind

declare i32 @get_global_size(i32) nounwind

declare i32 @get_local_size(i32) nounwind

declare i32 @get_group_id(i32) nounwind

define %struct.ungappedExtension addrspace(1)* @ungappedExtension_oneHitExtendD(i8* %sequenceStart, <4 x i32> %queryOffset, i8* %subjectHit, %struct.PSSMatrixFP* byval %PSSMatrixFP, i8* %subject, <4 x i32>* %sequenceHitEnd, i8 zeroext %encoding_numCodes, <4 x i32> %statistics_ungappedNominalDropoff, <4 x i32> %blast_ungappedNominalTrigger, %struct.ungappedExtension* %ungappedExtension_extensions, <4 x i32>* %numOfTriggerExtensions, <4 x i32> %sequenceCount, <4 x i32> %tid) nounwind {
entry:
  %retval = alloca %struct.ungappedExtension addrspace(1)*, align 4 ; <%struct.ungappedExtension addrspace(1)**> [#uses=3]
  %sequenceStart.addr = alloca i8*, align 4       ; <i8**> [#uses=2]
  %queryOffset.addr = alloca <4 x i32>, align 16  ; <<4 x i32>*> [#uses=4]
  %subjectHit.addr = alloca i8*, align 4          ; <i8**> [#uses=4]
  %subject.addr = alloca i8*, align 4             ; <i8**> [#uses=5]
  %sequenceHitEnd.addr = alloca <4 x i32>*, align 4 ; <<4 x i32>**> [#uses=2]
  %encoding_numCodes.addr = alloca i8, align 1    ; <i8*> [#uses=8]
  %statistics_ungappedNominalDropoff.addr = alloca <4 x i32>, align 16 ; <<4 x i32>*> [#uses=2]
  %blast_ungappedNominalTrigger.addr = alloca <4 x i32>, align 16 ; <<4 x i32>*> [#uses=2]
  %ungappedExtension_extensions.addr = alloca %struct.ungappedExtension*, align 4 ; <%struct.ungappedExtension**> [#uses=2]
  %numOfTriggerExtensions.addr = alloca <4 x i32>*, align 4 ; <<4 x i32>**> [#uses=3]
  %sequenceCount.addr = alloca <4 x i32>, align 16 ; <<4 x i32>*> [#uses=2]
  %tid.addr = alloca <4 x i32>, align 16          ; <<4 x i32>*> [#uses=1]
  %queryPosition = alloca <2 x i32>*, align 4     ; <<2 x i32>**> [#uses=14]
  %subjectPosition = alloca i8*, align 4          ; <i8**> [#uses=17]
  %subjectStart = alloca i8*, align 4             ; <i8**> [#uses=5]
  %subjectEnd = alloca i8*, align 4               ; <i8**> [#uses=5]
  %dropoff = alloca <4 x i32>, align 16           ; <<4 x i32>*> [#uses=6]
  %originalDropoff = alloca <4 x i32>, align 16   ; <<4 x i32>*> [#uses=4]
  %ungappedExtension_bestScore = alloca <4 x i32>, align 16 ; <<4 x i32>*> [#uses=10]
  %changeSinceBest = alloca <4 x i32>, align 16   ; <<4 x i32>*> [#uses=17]
  %tmp64 = alloca <4 x i32>, align 16             ; <<4 x i32>*> [#uses=2]
  %tmp143 = alloca <4 x i32>, align 1             ; <<4 x i32>*> [#uses=2]
  %diagonal = alloca <4 x i32>, align 16          ; <<4 x i32>*> [#uses=3]
  %newUngappedExtension = alloca %struct.ungappedExtension*, align 4 ; <%struct.ungappedExtension**> [#uses=15]
  %tmp234 = alloca %struct.coordinate, align 16   ; <%struct.coordinate*> [#uses=2]
  %tmp237 = alloca %struct.PSSMatrixFP, align 16  ; <%struct.PSSMatrixFP*> [#uses=2]
  store i8* %sequenceStart, i8** %sequenceStart.addr
  store <4 x i32> %queryOffset, <4 x i32>* %queryOffset.addr
  store i8* %subjectHit, i8** %subjectHit.addr
  store i8* %subject, i8** %subject.addr
  store <4 x i32>* %sequenceHitEnd, <4 x i32>** %sequenceHitEnd.addr
  store i8 %encoding_numCodes, i8* %encoding_numCodes.addr
  store <4 x i32> %statistics_ungappedNominalDropoff, <4 x i32>* %statistics_ungappedNominalDropoff.addr
  store <4 x i32> %blast_ungappedNominalTrigger, <4 x i32>* %blast_ungappedNominalTrigger.addr
  store %struct.ungappedExtension* %ungappedExtension_extensions, %struct.ungappedExtension** %ungappedExtension_extensions.addr
  store <4 x i32>* %numOfTriggerExtensions, <4 x i32>** %numOfTriggerExtensions.addr
  store <4 x i32> %sequenceCount, <4 x i32>* %sequenceCount.addr
  store <4 x i32> %tid, <4 x i32>* %tid.addr
  store <4 x i32> zeroinitializer, <4 x i32>* %changeSinceBest
  %tmp = load <4 x i32>* %statistics_ungappedNominalDropoff.addr ; <<4 x i32>> [#uses=1]
  %tmp1 = sub <4 x i32> zeroinitializer, %tmp     ; <<4 x i32>> [#uses=2]
  store <4 x i32> %tmp1, <4 x i32>* %dropoff
  store <4 x i32> %tmp1, <4 x i32>* %originalDropoff
  store <4 x i32> zeroinitializer, <4 x i32>* %ungappedExtension_bestScore
  %structele = getelementptr inbounds %struct.PSSMatrixFP* %PSSMatrixFP, i32 0, i32 4 ; <<2 x i32>**> [#uses=1]
  %tmp2 = load <2 x i32>** %structele             ; <<2 x i32>*> [#uses=1]
  %tmp3 = load <4 x i32>* %queryOffset.addr       ; <<4 x i32>> [#uses=1]
  %tmp4 = extractelement <4 x i32> %tmp3, i32 0   ; <i32> [#uses=1]
  %tmp5 = load i8* %encoding_numCodes.addr        ; <i8> [#uses=1]
  %conv = zext i8 %tmp5 to i32                    ; <i32> [#uses=1]
  %tmp6 = mul i32 %tmp4, %conv                    ; <i32> [#uses=1]
  %tmp7 = getelementptr <2 x i32>* %tmp2, i32 %tmp6 ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp7, <2 x i32>** %queryPosition
  %tmp8 = load i8** %subjectHit.addr              ; <i8*> [#uses=2]
  store i8* %tmp8, i8** %subjectStart
  store i8* %tmp8, i8** %subjectPosition
  br label %while.cond

return:                                           ; preds = %if.else, %if.then173
  %tmp257 = load %struct.ungappedExtension addrspace(1)** %retval ; <%struct.ungappedExtension addrspace(1)*> [#uses=1]
  ret %struct.ungappedExtension addrspace(1)* %tmp257

while.cond:                                       ; preds = %if.end, %entry
  %tmp9 = load <4 x i32>* %changeSinceBest        ; <<4 x i32>> [#uses=1]
  %tmp10 = extractelement <4 x i32> %tmp9, i32 0  ; <i32> [#uses=1]
  %tmp11 = load <4 x i32>* %dropoff               ; <<4 x i32>> [#uses=1]
  %tmp12 = extractelement <4 x i32> %tmp11, i32 0 ; <i32> [#uses=1]
  %cmp = icmp sgt i32 %tmp10, %tmp12              ; <i1> [#uses=1]
  br i1 %cmp, label %while.body, label %while.exit

while.exit:                                       ; preds = %while.cond
  %tmp58 = load i8** %subjectStart                ; <i8*> [#uses=1]
  %tmp59 = getelementptr i8* %tmp58, i32 1        ; <i8*> [#uses=1]
  store i8* %tmp59, i8** %subjectStart
  %structele60 = getelementptr inbounds %struct.PSSMatrixFP* %PSSMatrixFP, i32 0, i32 4 ; <<2 x i32>**> [#uses=1]
  %tmp61 = load <2 x i32>** %structele60          ; <<2 x i32>*> [#uses=1]
  %tmp62 = load <4 x i32>* %queryOffset.addr      ; <<4 x i32>> [#uses=1]
  %tmp63 = add <4 x i32> %tmp62, <i32 1, i32 1, i32 1, i32 1> ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp63, <4 x i32>* %tmp64
  %tmp65 = load <4 x i32>* %tmp64                 ; <<4 x i32>> [#uses=1]
  %tmp66 = extractelement <4 x i32> %tmp65, i32 0 ; <i32> [#uses=1]
  %tmp67 = load i8* %encoding_numCodes.addr       ; <i8> [#uses=1]
  %conv68 = zext i8 %tmp67 to i32                 ; <i32> [#uses=1]
  %tmp69 = mul i32 %tmp66, %conv68                ; <i32> [#uses=1]
  %tmp70 = getelementptr <2 x i32>* %tmp61, i32 %tmp69 ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp70, <2 x i32>** %queryPosition
  %tmp71 = load i8** %subjectHit.addr             ; <i8*> [#uses=1]
  %tmp72 = getelementptr i8* %tmp71, i32 1        ; <i8*> [#uses=2]
  store i8* %tmp72, i8** %subjectEnd
  store i8* %tmp72, i8** %subjectPosition
  store <4 x i32> zeroinitializer, <4 x i32>* %changeSinceBest
  %tmp73 = load <4 x i32>* %ungappedExtension_bestScore ; <<4 x i32>> [#uses=1]
  %tmp74 = extractelement <4 x i32> %tmp73, i32 0 ; <i32> [#uses=1]
  %tmp75 = sub i32 0, %tmp74                      ; <i32> [#uses=1]
  %tmp76 = load <4 x i32>* %originalDropoff       ; <<4 x i32>> [#uses=1]
  %tmp77 = extractelement <4 x i32> %tmp76, i32 0 ; <i32> [#uses=1]
  %cmp78 = icmp sgt i32 %tmp75, %tmp77            ; <i1> [#uses=1]
  br i1 %cmp78, label %if.then80, label %if.end79

while.body:                                       ; preds = %while.cond
  %tmp13 = load <4 x i32>* %changeSinceBest       ; <<4 x i32>> [#uses=1]
  %tmp14 = extractelement <4 x i32> %tmp13, i32 0 ; <i32> [#uses=1]
  %tmp15 = load <2 x i32>** %queryPosition        ; <<2 x i32>*> [#uses=1]
  %tmp16 = load i8** %subjectPosition             ; <i8*> [#uses=1]
  %tmp17 = load i8* %tmp16                        ; <i8> [#uses=1]
  %conv18 = zext i8 %tmp17 to i32                 ; <i32> [#uses=1]
  %arrayidx = getelementptr <2 x i32>* %tmp15, i32 %conv18 ; <<2 x i32>*> [#uses=1]
  %tmp19 = load <2 x i32>* %arrayidx              ; <<2 x i32>> [#uses=1]
  %tmp20 = extractelement <2 x i32> %tmp19, i32 0 ; <i32> [#uses=1]
  %tmp21 = add i32 %tmp14, %tmp20                 ; <i32> [#uses=1]
  %tmp22 = load <4 x i32>* %changeSinceBest       ; <<4 x i32>> [#uses=1]
  %tmp23 = insertelement <4 x i32> %tmp22, i32 %tmp21, i32 0 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp23, <4 x i32>* %changeSinceBest
  %tmp24 = load <4 x i32>* %changeSinceBest       ; <<4 x i32>> [#uses=1]
  %tmp25 = extractelement <4 x i32> %tmp24, i32 0 ; <i32> [#uses=1]
  %cmp26 = icmp sgt i32 %tmp25, 0                 ; <i1> [#uses=1]
  br i1 %cmp26, label %if.then, label %if.end

if.end:                                           ; preds = %do.exit, %while.body
  %tmp51 = load <2 x i32>** %queryPosition        ; <<2 x i32>*> [#uses=1]
  %tmp52 = load i8* %encoding_numCodes.addr       ; <i8> [#uses=1]
  %conv53 = zext i8 %tmp52 to i32                 ; <i32> [#uses=1]
  %sub.ptr.neg54 = sub i32 0, %conv53             ; <i32> [#uses=1]
  %tmp55 = getelementptr <2 x i32>* %tmp51, i32 %sub.ptr.neg54 ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp55, <2 x i32>** %queryPosition
  %tmp56 = load i8** %subjectPosition             ; <i8*> [#uses=1]
  %tmp57 = getelementptr i8* %tmp56, i32 -1       ; <i8*> [#uses=1]
  store i8* %tmp57, i8** %subjectPosition
  br label %while.cond

if.then:                                          ; preds = %while.body
  br label %do.body

do.exit:                                          ; preds = %do.cond
  %tmp50 = load i8** %subjectPosition             ; <i8*> [#uses=1]
  store i8* %tmp50, i8** %subjectStart
  br label %if.end

do.body:                                          ; preds = %do.cond, %if.then
  %tmp27 = load <4 x i32>* %ungappedExtension_bestScore ; <<4 x i32>> [#uses=1]
  %tmp28 = load <4 x i32>* %changeSinceBest       ; <<4 x i32>> [#uses=1]
  %tmp29 = add <4 x i32> %tmp27, %tmp28           ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp29, <4 x i32>* %ungappedExtension_bestScore
  %tmp30 = load <2 x i32>** %queryPosition        ; <<2 x i32>*> [#uses=1]
  %tmp31 = load i8* %encoding_numCodes.addr       ; <i8> [#uses=1]
  %conv32 = zext i8 %tmp31 to i32                 ; <i32> [#uses=1]
  %sub.ptr.neg = sub i32 0, %conv32               ; <i32> [#uses=1]
  %tmp33 = getelementptr <2 x i32>* %tmp30, i32 %sub.ptr.neg ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp33, <2 x i32>** %queryPosition
  %tmp34 = load i8** %subjectPosition             ; <i8*> [#uses=1]
  %tmp35 = getelementptr i8* %tmp34, i32 -1       ; <i8*> [#uses=1]
  store i8* %tmp35, i8** %subjectPosition
  %tmp36 = load <2 x i32>** %queryPosition        ; <<2 x i32>*> [#uses=1]
  %tmp37 = load i8** %subjectPosition             ; <i8*> [#uses=1]
  %tmp38 = load i8* %tmp37                        ; <i8> [#uses=1]
  %conv39 = zext i8 %tmp38 to i32                 ; <i32> [#uses=1]
  %arrayidx40 = getelementptr <2 x i32>* %tmp36, i32 %conv39 ; <<2 x i32>*> [#uses=1]
  %tmp41 = load <2 x i32>* %arrayidx40            ; <<2 x i32>> [#uses=1]
  %tmp42 = extractelement <2 x i32> %tmp41, i32 0 ; <i32> [#uses=4]
  %conv43 = insertelement <4 x i32> undef, i32 %tmp42, i32 0 ; <<4 x i32>> [#uses=1]
  %conv44 = insertelement <4 x i32> %conv43, i32 %tmp42, i32 1 ; <<4 x i32>> [#uses=1]
  %conv45 = insertelement <4 x i32> %conv44, i32 %tmp42, i32 2 ; <<4 x i32>> [#uses=1]
  %conv46 = insertelement <4 x i32> %conv45, i32 %tmp42, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv46, <4 x i32>* %changeSinceBest
  br label %do.cond

do.cond:                                          ; preds = %do.body
  %tmp47 = load <4 x i32>* %changeSinceBest       ; <<4 x i32>> [#uses=1]
  %tmp48 = extractelement <4 x i32> %tmp47, i32 0 ; <i32> [#uses=1]
  %cmp49 = icmp sgt i32 %tmp48, 0                 ; <i1> [#uses=1]
  br i1 %cmp49, label %do.body, label %do.exit

if.end79:                                         ; preds = %if.then80, %while.exit
  br label %while.cond83

if.then80:                                        ; preds = %while.exit
  %tmp81 = load <4 x i32>* %ungappedExtension_bestScore ; <<4 x i32>> [#uses=1]
  %tmp82 = sub <4 x i32> zeroinitializer, %tmp81  ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp82, <4 x i32>* %dropoff
  br label %if.end79

while.cond83:                                     ; preds = %if.end107, %if.end79
  %tmp84 = load <4 x i32>* %changeSinceBest       ; <<4 x i32>> [#uses=1]
  %tmp85 = extractelement <4 x i32> %tmp84, i32 0 ; <i32> [#uses=1]
  %tmp86 = load <4 x i32>* %dropoff               ; <<4 x i32>> [#uses=1]
  %tmp87 = extractelement <4 x i32> %tmp86, i32 0 ; <i32> [#uses=1]
  %cmp88 = icmp sgt i32 %tmp85, %tmp87            ; <i1> [#uses=1]
  br i1 %cmp88, label %while.body90, label %while.exit89

while.exit89:                                     ; preds = %while.cond83
  %tmp158 = load i8** %subjectEnd                 ; <i8*> [#uses=1]
  %tmp159 = getelementptr i8* %tmp158, i32 -1     ; <i8*> [#uses=1]
  store i8* %tmp159, i8** %subjectEnd
  %tmp160 = load <4 x i32>** %sequenceHitEnd.addr ; <<4 x i32>*> [#uses=1]
  %tmp161 = load i8** %subjectPosition            ; <i8*> [#uses=1]
  %tmp162 = load i8** %sequenceStart.addr         ; <i8*> [#uses=1]
  %sub.ptr.lhs.cast = ptrtoint i8* %tmp161 to i32 ; <i32> [#uses=1]
  %sub.ptr.rhs.cast = ptrtoint i8* %tmp162 to i32 ; <i32> [#uses=1]
  %sub.ptr.sub = sub i32 %sub.ptr.lhs.cast, %sub.ptr.rhs.cast ; <i32> [#uses=4]
  %conv163 = insertelement <4 x i32> undef, i32 %sub.ptr.sub, i32 0 ; <<4 x i32>> [#uses=1]
  %conv164 = insertelement <4 x i32> %conv163, i32 %sub.ptr.sub, i32 1 ; <<4 x i32>> [#uses=1]
  %conv165 = insertelement <4 x i32> %conv164, i32 %sub.ptr.sub, i32 2 ; <<4 x i32>> [#uses=1]
  %conv166 = insertelement <4 x i32> %conv165, i32 %sub.ptr.sub, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv166, <4 x i32>* %tmp160
  %tmp167 = load <4 x i32>* %ungappedExtension_bestScore ; <<4 x i32>> [#uses=1]
  %tmp168 = extractelement <4 x i32> %tmp167, i32 0 ; <i32> [#uses=1]
  %tmp169 = load <4 x i32>* %blast_ungappedNominalTrigger.addr ; <<4 x i32>> [#uses=1]
  %tmp170 = extractelement <4 x i32> %tmp169, i32 0 ; <i32> [#uses=1]
  %cmp171 = icmp sge i32 %tmp168, %tmp170         ; <i1> [#uses=1]
  br i1 %cmp171, label %if.then173, label %if.else

while.body90:                                     ; preds = %while.cond83
  %tmp91 = load <4 x i32>* %changeSinceBest       ; <<4 x i32>> [#uses=1]
  %tmp92 = load <2 x i32>** %queryPosition        ; <<2 x i32>*> [#uses=1]
  %tmp93 = load i8** %subjectPosition             ; <i8*> [#uses=1]
  %tmp94 = load i8* %tmp93                        ; <i8> [#uses=1]
  %conv95 = zext i8 %tmp94 to i32                 ; <i32> [#uses=1]
  %arrayidx96 = getelementptr <2 x i32>* %tmp92, i32 %conv95 ; <<2 x i32>*> [#uses=1]
  %tmp97 = load <2 x i32>* %arrayidx96            ; <<2 x i32>> [#uses=1]
  %tmp98 = extractelement <2 x i32> %tmp97, i32 0 ; <i32> [#uses=4]
  %conv99 = insertelement <4 x i32> undef, i32 %tmp98, i32 0 ; <<4 x i32>> [#uses=1]
  %conv100 = insertelement <4 x i32> %conv99, i32 %tmp98, i32 1 ; <<4 x i32>> [#uses=1]
  %conv101 = insertelement <4 x i32> %conv100, i32 %tmp98, i32 2 ; <<4 x i32>> [#uses=1]
  %conv102 = insertelement <4 x i32> %conv101, i32 %tmp98, i32 3 ; <<4 x i32>> [#uses=1]
  %tmp103 = add <4 x i32> %tmp91, %conv102        ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp103, <4 x i32>* %changeSinceBest
  %tmp104 = load <4 x i32>* %changeSinceBest      ; <<4 x i32>> [#uses=1]
  %tmp105 = extractelement <4 x i32> %tmp104, i32 0 ; <i32> [#uses=1]
  %cmp106 = icmp sgt i32 %tmp105, 0               ; <i1> [#uses=1]
  br i1 %cmp106, label %if.then108, label %if.end107

if.end107:                                        ; preds = %if.end149, %while.body90
  %tmp152 = load <2 x i32>** %queryPosition       ; <<2 x i32>*> [#uses=1]
  %tmp153 = load i8* %encoding_numCodes.addr      ; <i8> [#uses=1]
  %conv154 = zext i8 %tmp153 to i32               ; <i32> [#uses=1]
  %tmp155 = getelementptr <2 x i32>* %tmp152, i32 %conv154 ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp155, <2 x i32>** %queryPosition
  %tmp156 = load i8** %subjectPosition            ; <i8*> [#uses=1]
  %tmp157 = getelementptr i8* %tmp156, i32 1      ; <i8*> [#uses=1]
  store i8* %tmp157, i8** %subjectPosition
  br label %while.cond83

if.then108:                                       ; preds = %while.body90
  br label %do.body110

do.exit109:                                       ; preds = %do.cond111
  %tmp135 = load i8** %subjectPosition            ; <i8*> [#uses=1]
  store i8* %tmp135, i8** %subjectEnd
  %tmp136 = load <4 x i32>* %ungappedExtension_bestScore ; <<4 x i32>> [#uses=1]
  %tmp137 = extractelement <4 x i32> %tmp136, i32 0 ; <i32> [#uses=1]
  %tmp138 = sub i32 0, %tmp137                    ; <i32> [#uses=4]
  %conv139 = insertelement <4 x i32> undef, i32 %tmp138, i32 0 ; <<4 x i32>> [#uses=1]
  %conv140 = insertelement <4 x i32> %conv139, i32 %tmp138, i32 1 ; <<4 x i32>> [#uses=1]
  %conv141 = insertelement <4 x i32> %conv140, i32 %tmp138, i32 2 ; <<4 x i32>> [#uses=1]
  %conv142 = insertelement <4 x i32> %conv141, i32 %tmp138, i32 3 ; <<4 x i32>> [#uses=2]
  store <4 x i32> %conv142, <4 x i32>* %dropoff
  store <4 x i32> %conv142, <4 x i32>* %tmp143
  %tmp144 = load <4 x i32>* %tmp143               ; <<4 x i32>> [#uses=1]
  %tmp145 = extractelement <4 x i32> %tmp144, i32 0 ; <i32> [#uses=1]
  %tmp146 = load <4 x i32>* %originalDropoff      ; <<4 x i32>> [#uses=1]
  %tmp147 = extractelement <4 x i32> %tmp146, i32 0 ; <i32> [#uses=1]
  %cmp148 = icmp slt i32 %tmp145, %tmp147         ; <i1> [#uses=1]
  br i1 %cmp148, label %if.then150, label %if.end149

do.body110:                                       ; preds = %do.cond111, %if.then108
  %tmp112 = load <4 x i32>* %ungappedExtension_bestScore ; <<4 x i32>> [#uses=1]
  %tmp113 = load <4 x i32>* %changeSinceBest      ; <<4 x i32>> [#uses=1]
  %tmp114 = add <4 x i32> %tmp112, %tmp113        ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp114, <4 x i32>* %ungappedExtension_bestScore
  %tmp115 = load <2 x i32>** %queryPosition       ; <<2 x i32>*> [#uses=1]
  %tmp116 = load i8* %encoding_numCodes.addr      ; <i8> [#uses=1]
  %conv117 = zext i8 %tmp116 to i32               ; <i32> [#uses=1]
  %tmp118 = getelementptr <2 x i32>* %tmp115, i32 %conv117 ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp118, <2 x i32>** %queryPosition
  %tmp119 = load i8** %subjectPosition            ; <i8*> [#uses=1]
  %tmp120 = getelementptr i8* %tmp119, i32 1      ; <i8*> [#uses=1]
  store i8* %tmp120, i8** %subjectPosition
  %tmp121 = load <2 x i32>** %queryPosition       ; <<2 x i32>*> [#uses=1]
  %tmp122 = load i8** %subjectPosition            ; <i8*> [#uses=1]
  %tmp123 = load i8* %tmp122                      ; <i8> [#uses=1]
  %conv124 = zext i8 %tmp123 to i32               ; <i32> [#uses=1]
  %arrayidx125 = getelementptr <2 x i32>* %tmp121, i32 %conv124 ; <<2 x i32>*> [#uses=1]
  %tmp126 = load <2 x i32>* %arrayidx125          ; <<2 x i32>> [#uses=1]
  %tmp127 = extractelement <2 x i32> %tmp126, i32 0 ; <i32> [#uses=4]
  %conv128 = insertelement <4 x i32> undef, i32 %tmp127, i32 0 ; <<4 x i32>> [#uses=1]
  %conv129 = insertelement <4 x i32> %conv128, i32 %tmp127, i32 1 ; <<4 x i32>> [#uses=1]
  %conv130 = insertelement <4 x i32> %conv129, i32 %tmp127, i32 2 ; <<4 x i32>> [#uses=1]
  %conv131 = insertelement <4 x i32> %conv130, i32 %tmp127, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv131, <4 x i32>* %changeSinceBest
  br label %do.cond111

do.cond111:                                       ; preds = %do.body110
  %tmp132 = load <4 x i32>* %changeSinceBest      ; <<4 x i32>> [#uses=1]
  %tmp133 = extractelement <4 x i32> %tmp132, i32 0 ; <i32> [#uses=1]
  %cmp134 = icmp sgt i32 %tmp133, 0               ; <i1> [#uses=1]
  br i1 %cmp134, label %do.body110, label %do.exit109

if.end149:                                        ; preds = %if.then150, %do.exit109
  br label %if.end107

if.then150:                                       ; preds = %do.exit109
  %tmp151 = load <4 x i32>* %originalDropoff      ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp151, <4 x i32>* %dropoff
  br label %if.end149

if.then173:                                       ; preds = %while.exit89
  store %struct.ungappedExtension* null, %struct.ungappedExtension** %newUngappedExtension
  %tmp174 = load %struct.ungappedExtension** %ungappedExtension_extensions.addr ; <%struct.ungappedExtension*> [#uses=1]
  %tmp175 = load <4 x i32>** %numOfTriggerExtensions.addr ; <<4 x i32>*> [#uses=1]
  %tmp176 = load <4 x i32>* %tmp175               ; <<4 x i32>> [#uses=1]
  %tmp177 = extractelement <4 x i32> %tmp176, i32 0 ; <i32> [#uses=1]
  %tmp178 = getelementptr %struct.ungappedExtension* %tmp174, i32 %tmp177 ; <%struct.ungappedExtension*> [#uses=1]
  store %struct.ungappedExtension* %tmp178, %struct.ungappedExtension** %newUngappedExtension
  %tmp179 = load i8** %subjectHit.addr            ; <i8*> [#uses=1]
  %tmp180 = load i8** %subject.addr               ; <i8*> [#uses=1]
  %sub.ptr.lhs.cast181 = ptrtoint i8* %tmp179 to i32 ; <i32> [#uses=1]
  %sub.ptr.rhs.cast182 = ptrtoint i8* %tmp180 to i32 ; <i32> [#uses=1]
  %sub.ptr.sub183 = sub i32 %sub.ptr.lhs.cast181, %sub.ptr.rhs.cast182 ; <i32> [#uses=4]
  %conv184 = insertelement <4 x i32> undef, i32 %sub.ptr.sub183, i32 0 ; <<4 x i32>> [#uses=1]
  %conv185 = insertelement <4 x i32> %conv184, i32 %sub.ptr.sub183, i32 1 ; <<4 x i32>> [#uses=1]
  %conv186 = insertelement <4 x i32> %conv185, i32 %sub.ptr.sub183, i32 2 ; <<4 x i32>> [#uses=1]
  %conv187 = insertelement <4 x i32> %conv186, i32 %sub.ptr.sub183, i32 3 ; <<4 x i32>> [#uses=1]
  %tmp188 = load <4 x i32>* %queryOffset.addr     ; <<4 x i32>> [#uses=1]
  %tmp189 = sub <4 x i32> %conv187, %tmp188       ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp189, <4 x i32>* %diagonal
  %tmp190 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele191 = getelementptr inbounds %struct.ungappedExtension* %tmp190, i32 0, i32 0 ; <%struct.coordinate*> [#uses=1]
  %structele192 = getelementptr inbounds %struct.coordinate* %structele191, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp193 = load i8** %subjectStart               ; <i8*> [#uses=1]
  %tmp194 = load i8** %subject.addr               ; <i8*> [#uses=1]
  %sub.ptr.lhs.cast195 = ptrtoint i8* %tmp193 to i32 ; <i32> [#uses=1]
  %sub.ptr.rhs.cast196 = ptrtoint i8* %tmp194 to i32 ; <i32> [#uses=1]
  %sub.ptr.sub197 = sub i32 %sub.ptr.lhs.cast195, %sub.ptr.rhs.cast196 ; <i32> [#uses=4]
  %conv198 = insertelement <4 x i32> undef, i32 %sub.ptr.sub197, i32 0 ; <<4 x i32>> [#uses=1]
  %conv199 = insertelement <4 x i32> %conv198, i32 %sub.ptr.sub197, i32 1 ; <<4 x i32>> [#uses=1]
  %conv200 = insertelement <4 x i32> %conv199, i32 %sub.ptr.sub197, i32 2 ; <<4 x i32>> [#uses=1]
  %conv201 = insertelement <4 x i32> %conv200, i32 %sub.ptr.sub197, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv201, <4 x i32>* %structele192
  %tmp202 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele203 = getelementptr inbounds %struct.ungappedExtension* %tmp202, i32 0, i32 1 ; <%struct.coordinate*> [#uses=1]
  %structele204 = getelementptr inbounds %struct.coordinate* %structele203, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp205 = load i8** %subjectEnd                 ; <i8*> [#uses=1]
  %tmp206 = load i8** %subject.addr               ; <i8*> [#uses=1]
  %sub.ptr.lhs.cast207 = ptrtoint i8* %tmp205 to i32 ; <i32> [#uses=1]
  %sub.ptr.rhs.cast208 = ptrtoint i8* %tmp206 to i32 ; <i32> [#uses=1]
  %sub.ptr.sub209 = sub i32 %sub.ptr.lhs.cast207, %sub.ptr.rhs.cast208 ; <i32> [#uses=4]
  %conv210 = insertelement <4 x i32> undef, i32 %sub.ptr.sub209, i32 0 ; <<4 x i32>> [#uses=1]
  %conv211 = insertelement <4 x i32> %conv210, i32 %sub.ptr.sub209, i32 1 ; <<4 x i32>> [#uses=1]
  %conv212 = insertelement <4 x i32> %conv211, i32 %sub.ptr.sub209, i32 2 ; <<4 x i32>> [#uses=1]
  %conv213 = insertelement <4 x i32> %conv212, i32 %sub.ptr.sub209, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv213, <4 x i32>* %structele204
  %tmp214 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele215 = getelementptr inbounds %struct.ungappedExtension* %tmp214, i32 0, i32 0 ; <%struct.coordinate*> [#uses=1]
  %structele216 = getelementptr inbounds %struct.coordinate* %structele215, i32 0, i32 0 ; <<4 x i32>*> [#uses=1]
  %tmp217 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele218 = getelementptr inbounds %struct.ungappedExtension* %tmp217, i32 0, i32 0 ; <%struct.coordinate*> [#uses=1]
  %structele219 = getelementptr inbounds %struct.coordinate* %structele218, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp220 = load <4 x i32>* %structele219         ; <<4 x i32>> [#uses=1]
  %tmp221 = load <4 x i32>* %diagonal             ; <<4 x i32>> [#uses=1]
  %tmp222 = sub <4 x i32> %tmp220, %tmp221        ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp222, <4 x i32>* %structele216
  %tmp223 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele224 = getelementptr inbounds %struct.ungappedExtension* %tmp223, i32 0, i32 1 ; <%struct.coordinate*> [#uses=1]
  %structele225 = getelementptr inbounds %struct.coordinate* %structele224, i32 0, i32 0 ; <<4 x i32>*> [#uses=1]
  %tmp226 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele227 = getelementptr inbounds %struct.ungappedExtension* %tmp226, i32 0, i32 1 ; <%struct.coordinate*> [#uses=1]
  %structele228 = getelementptr inbounds %struct.coordinate* %structele227, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp229 = load <4 x i32>* %structele228         ; <<4 x i32>> [#uses=1]
  %tmp230 = load <4 x i32>* %diagonal             ; <<4 x i32>> [#uses=1]
  %tmp231 = sub <4 x i32> %tmp229, %tmp230        ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp231, <4 x i32>* %structele225
  %tmp232 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele233 = getelementptr inbounds %struct.ungappedExtension* %tmp232, i32 0, i32 2 ; <%struct.coordinate*> [#uses=1]
  %tmp235 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %tmp236 = load %struct.PSSMatrixFP* %PSSMatrixFP ; <%struct.PSSMatrixFP> [#uses=1]
  store %struct.PSSMatrixFP %tmp236, %struct.PSSMatrixFP* %tmp237
  %tmp238 = load i8** %subject.addr               ; <i8*> [#uses=1]
  %tmp239 = load i8* %encoding_numCodes.addr      ; <i8> [#uses=1]
  call void @ungappedExtension_findProteinSeed(%struct.coordinate* noalias sret %tmp234, %struct.ungappedExtension* %tmp235, %struct.PSSMatrixFP* byval %tmp237, i8* %tmp238, i8 zeroext %tmp239) nounwind
  %tmp240 = load %struct.coordinate* %tmp234      ; <%struct.coordinate> [#uses=1]
  store %struct.coordinate %tmp240, %struct.coordinate* %structele233
  %tmp241 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele242 = getelementptr inbounds %struct.ungappedExtension* %tmp241, i32 0, i32 10 ; <%struct.ungappedExtension**> [#uses=1]
  store %struct.ungappedExtension* null, %struct.ungappedExtension** %structele242
  %tmp243 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele244 = getelementptr inbounds %struct.ungappedExtension* %tmp243, i32 0, i32 3 ; <<4 x i32>*> [#uses=1]
  %tmp245 = load <4 x i32>* %ungappedExtension_bestScore ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp245, <4 x i32>* %structele244
  %tmp246 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele247 = getelementptr inbounds %struct.ungappedExtension* %tmp246, i32 0, i32 6 ; <i8*> [#uses=1]
  store i8 1, i8* %structele247
  %tmp248 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele249 = getelementptr inbounds %struct.ungappedExtension* %tmp248, i32 0, i32 4 ; <<4 x i32>*> [#uses=1]
  %tmp250 = load <4 x i32>* %sequenceCount.addr   ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp250, <4 x i32>* %structele249
  %tmp251 = load <4 x i32>** %numOfTriggerExtensions.addr ; <<4 x i32>*> [#uses=2]
  %tmp252 = load <4 x i32>* %tmp251               ; <<4 x i32>> [#uses=1]
  %tmp253 = add <4 x i32> %tmp252, <i32 1, i32 1, i32 1, i32 1> ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp253, <4 x i32>* %tmp251
  %tmp254 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %conv255 = ptrtoint %struct.ungappedExtension* %tmp254 to i32 ; <i32> [#uses=1]
  %conv256 = inttoptr i32 %conv255 to %struct.ungappedExtension addrspace(1)* ; <%struct.ungappedExtension addrspace(1)*> [#uses=1]
  store %struct.ungappedExtension addrspace(1)* %conv256, %struct.ungappedExtension addrspace(1)** %retval
  br label %return

if.else:                                          ; preds = %while.exit89
  store %struct.ungappedExtension addrspace(1)* null, %struct.ungappedExtension addrspace(1)** %retval
  br label %return
}

define void @ungappedExtension_findProteinSeed(%struct.coordinate* noalias sret %agg.result, %struct.ungappedExtension* %ungappedExtension, %struct.PSSMatrixFP* byval %PSSMatrixFP, i8* %subject, i8 zeroext %encoding_numCodes) nounwind {
entry:
  %retval = alloca %struct.coordinate, align 16   ; <%struct.coordinate*> [#uses=2]
  %ungappedExtension.addr = alloca %struct.ungappedExtension*, align 4 ; <%struct.ungappedExtension**> [#uses=10]
  %subject.addr = alloca i8*, align 4             ; <i8**> [#uses=3]
  %encoding_numCodes.addr = alloca i8, align 1    ; <i8*> [#uses=8]
  %queryWindowStart = alloca <2 x i32>*, align 4  ; <<2 x i32>**> [#uses=6]
  %queryWindowEnd = alloca <2 x i32>*, align 4    ; <<2 x i32>**> [#uses=10]
  %subjectWindowStart = alloca i8*, align 4       ; <i8**> [#uses=6]
  %subjectWindowEnd = alloca i8*, align 4         ; <i8**> [#uses=9]
  %bestQueryPosition = alloca <2 x i32>*, align 4 ; <<2 x i32>**> [#uses=3]
  %bestSubjectPosition = alloca i8*, align 4      ; <i8**> [#uses=3]
  %bestSegmentScore = alloca <4 x i32>, align 16  ; <<4 x i32>*> [#uses=3]
  %nominalScore = alloca <4 x i32>, align 16      ; <<4 x i32>*> [#uses=10]
  %count = alloca <4 x i32>, align 16             ; <<4 x i32>*> [#uses=4]
  %seed = alloca %struct.coordinate, align 16     ; <%struct.coordinate*> [#uses=5]
  store %struct.ungappedExtension* %ungappedExtension, %struct.ungappedExtension** %ungappedExtension.addr
  store i8* %subject, i8** %subject.addr
  store i8 %encoding_numCodes, i8* %encoding_numCodes.addr
  %tmp = load %struct.ungappedExtension** %ungappedExtension.addr ; <%struct.ungappedExtension*> [#uses=1]
  %structele = getelementptr inbounds %struct.ungappedExtension* %tmp, i32 0, i32 1 ; <%struct.coordinate*> [#uses=1]
  %structele1 = getelementptr inbounds %struct.coordinate* %structele, i32 0, i32 0 ; <<4 x i32>*> [#uses=1]
  %tmp2 = load <4 x i32>* %structele1             ; <<4 x i32>> [#uses=1]
  %tmp3 = extractelement <4 x i32> %tmp2, i32 0   ; <i32> [#uses=1]
  %tmp4 = load %struct.ungappedExtension** %ungappedExtension.addr ; <%struct.ungappedExtension*> [#uses=1]
  %structele5 = getelementptr inbounds %struct.ungappedExtension* %tmp4, i32 0, i32 0 ; <%struct.coordinate*> [#uses=1]
  %structele6 = getelementptr inbounds %struct.coordinate* %structele5, i32 0, i32 0 ; <<4 x i32>*> [#uses=1]
  %tmp7 = load <4 x i32>* %structele6             ; <<4 x i32>> [#uses=1]
  %tmp8 = extractelement <4 x i32> %tmp7, i32 0   ; <i32> [#uses=1]
  %tmp9 = sub i32 %tmp3, %tmp8                    ; <i32> [#uses=1]
  %cmp = icmp slt i32 %tmp9, 11                   ; <i1> [#uses=1]
  br i1 %cmp, label %if.then, label %if.else

return:                                           ; preds = %if.end
  %tmp169 = load %struct.coordinate* %retval      ; <%struct.coordinate> [#uses=1]
  store %struct.coordinate %tmp169, %struct.coordinate* %agg.result
  ret void

if.end:                                           ; preds = %while.exit95, %if.then
  %tmp168 = load %struct.coordinate* %seed        ; <%struct.coordinate> [#uses=1]
  store %struct.coordinate %tmp168, %struct.coordinate* %retval
  br label %return

if.then:                                          ; preds = %entry
  %structele10 = getelementptr inbounds %struct.coordinate* %seed, i32 0, i32 0 ; <<4 x i32>*> [#uses=1]
  %tmp11 = load %struct.ungappedExtension** %ungappedExtension.addr ; <%struct.ungappedExtension*> [#uses=1]
  %structele12 = getelementptr inbounds %struct.ungappedExtension* %tmp11, i32 0, i32 1 ; <%struct.coordinate*> [#uses=1]
  %structele13 = getelementptr inbounds %struct.coordinate* %structele12, i32 0, i32 0 ; <<4 x i32>*> [#uses=1]
  %tmp14 = load <4 x i32>* %structele13           ; <<4 x i32>> [#uses=1]
  %tmp15 = load %struct.ungappedExtension** %ungappedExtension.addr ; <%struct.ungappedExtension*> [#uses=1]
  %structele16 = getelementptr inbounds %struct.ungappedExtension* %tmp15, i32 0, i32 0 ; <%struct.coordinate*> [#uses=1]
  %structele17 = getelementptr inbounds %struct.coordinate* %structele16, i32 0, i32 0 ; <<4 x i32>*> [#uses=1]
  %tmp18 = load <4 x i32>* %structele17           ; <<4 x i32>> [#uses=1]
  %tmp19 = add <4 x i32> %tmp14, %tmp18           ; <<4 x i32>> [#uses=1]
  %tmp20 = sdiv <4 x i32> %tmp19, <i32 2, i32 2, i32 2, i32 2> ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp20, <4 x i32>* %structele10
  %structele21 = getelementptr inbounds %struct.coordinate* %seed, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp22 = load %struct.ungappedExtension** %ungappedExtension.addr ; <%struct.ungappedExtension*> [#uses=1]
  %structele23 = getelementptr inbounds %struct.ungappedExtension* %tmp22, i32 0, i32 1 ; <%struct.coordinate*> [#uses=1]
  %structele24 = getelementptr inbounds %struct.coordinate* %structele23, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp25 = load <4 x i32>* %structele24           ; <<4 x i32>> [#uses=1]
  %tmp26 = load %struct.ungappedExtension** %ungappedExtension.addr ; <%struct.ungappedExtension*> [#uses=1]
  %structele27 = getelementptr inbounds %struct.ungappedExtension* %tmp26, i32 0, i32 0 ; <%struct.coordinate*> [#uses=1]
  %structele28 = getelementptr inbounds %struct.coordinate* %structele27, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp29 = load <4 x i32>* %structele28           ; <<4 x i32>> [#uses=1]
  %tmp30 = add <4 x i32> %tmp25, %tmp29           ; <<4 x i32>> [#uses=1]
  %tmp31 = sdiv <4 x i32> %tmp30, <i32 2, i32 2, i32 2, i32 2> ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp31, <4 x i32>* %structele21
  br label %if.end

if.else:                                          ; preds = %entry
  %structele32 = getelementptr inbounds %struct.PSSMatrixFP* %PSSMatrixFP, i32 0, i32 4 ; <<2 x i32>**> [#uses=1]
  %tmp33 = load <2 x i32>** %structele32          ; <<2 x i32>*> [#uses=1]
  %tmp34 = load %struct.ungappedExtension** %ungappedExtension.addr ; <%struct.ungappedExtension*> [#uses=1]
  %structele35 = getelementptr inbounds %struct.ungappedExtension* %tmp34, i32 0, i32 0 ; <%struct.coordinate*> [#uses=1]
  %structele36 = getelementptr inbounds %struct.coordinate* %structele35, i32 0, i32 0 ; <<4 x i32>*> [#uses=1]
  %tmp37 = load <4 x i32>* %structele36           ; <<4 x i32>> [#uses=1]
  %tmp38 = extractelement <4 x i32> %tmp37, i32 0 ; <i32> [#uses=1]
  %tmp39 = load i8* %encoding_numCodes.addr       ; <i8> [#uses=1]
  %conv = zext i8 %tmp39 to i32                   ; <i32> [#uses=1]
  %tmp40 = mul i32 %tmp38, %conv                  ; <i32> [#uses=1]
  %tmp41 = getelementptr <2 x i32>* %tmp33, i32 %tmp40 ; <<2 x i32>*> [#uses=2]
  store <2 x i32>* %tmp41, <2 x i32>** %queryWindowEnd
  store <2 x i32>* %tmp41, <2 x i32>** %queryWindowStart
  %tmp42 = load i8** %subject.addr                ; <i8*> [#uses=1]
  %tmp43 = load %struct.ungappedExtension** %ungappedExtension.addr ; <%struct.ungappedExtension*> [#uses=1]
  %structele44 = getelementptr inbounds %struct.ungappedExtension* %tmp43, i32 0, i32 0 ; <%struct.coordinate*> [#uses=1]
  %structele45 = getelementptr inbounds %struct.coordinate* %structele44, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp46 = load <4 x i32>* %structele45           ; <<4 x i32>> [#uses=1]
  %tmp47 = extractelement <4 x i32> %tmp46, i32 0 ; <i32> [#uses=1]
  %tmp48 = getelementptr i8* %tmp42, i32 %tmp47   ; <i8*> [#uses=2]
  store i8* %tmp48, i8** %subjectWindowEnd
  store i8* %tmp48, i8** %subjectWindowStart
  store <4 x i32> zeroinitializer, <4 x i32>* %nominalScore
  store <4 x i32> zeroinitializer, <4 x i32>* %count
  br label %while.cond

while.cond:                                       ; preds = %while.body, %if.else
  %tmp49 = load <4 x i32>* %count                 ; <<4 x i32>> [#uses=1]
  %tmp50 = extractelement <4 x i32> %tmp49, i32 0 ; <i32> [#uses=1]
  %cmp51 = icmp slt i32 %tmp50, 11                ; <i1> [#uses=1]
  br i1 %cmp51, label %while.body, label %while.exit

while.exit:                                       ; preds = %while.cond
  %tmp72 = load i8* %encoding_numCodes.addr       ; <i8> [#uses=1]
  %tmp73 = load <2 x i32>** %queryWindowEnd       ; <<2 x i32>*> [#uses=1]
  %conv74 = zext i8 %tmp72 to i32                 ; <i32> [#uses=1]
  %sub.ptr.neg = sub i32 0, %conv74               ; <i32> [#uses=1]
  %tmp75 = getelementptr <2 x i32>* %tmp73, i32 %sub.ptr.neg ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp75, <2 x i32>** %queryWindowEnd
  %tmp76 = load i8** %subjectWindowEnd            ; <i8*> [#uses=1]
  %tmp77 = getelementptr i8* %tmp76, i32 -1       ; <i8*> [#uses=1]
  store i8* %tmp77, i8** %subjectWindowEnd
  %tmp78 = load <2 x i32>** %queryWindowStart     ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp78, <2 x i32>** %bestQueryPosition
  %tmp79 = load i8** %subjectWindowStart          ; <i8*> [#uses=1]
  store i8* %tmp79, i8** %bestSubjectPosition
  %tmp80 = load <4 x i32>* %nominalScore          ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp80, <4 x i32>* %bestSegmentScore
  br label %while.cond81

while.body:                                       ; preds = %while.cond
  %tmp52 = load <4 x i32>* %nominalScore          ; <<4 x i32>> [#uses=1]
  %tmp53 = load <2 x i32>** %queryWindowEnd       ; <<2 x i32>*> [#uses=1]
  %tmp54 = load i8** %subjectWindowEnd            ; <i8*> [#uses=1]
  %tmp55 = load i8* %tmp54                        ; <i8> [#uses=1]
  %conv56 = zext i8 %tmp55 to i32                 ; <i32> [#uses=1]
  %arrayidx = getelementptr <2 x i32>* %tmp53, i32 %conv56 ; <<2 x i32>*> [#uses=1]
  %tmp57 = load <2 x i32>* %arrayidx              ; <<2 x i32>> [#uses=1]
  %tmp58 = extractelement <2 x i32> %tmp57, i32 0 ; <i32> [#uses=4]
  %conv59 = insertelement <4 x i32> undef, i32 %tmp58, i32 0 ; <<4 x i32>> [#uses=1]
  %conv60 = insertelement <4 x i32> %conv59, i32 %tmp58, i32 1 ; <<4 x i32>> [#uses=1]
  %conv61 = insertelement <4 x i32> %conv60, i32 %tmp58, i32 2 ; <<4 x i32>> [#uses=1]
  %conv62 = insertelement <4 x i32> %conv61, i32 %tmp58, i32 3 ; <<4 x i32>> [#uses=1]
  %tmp63 = add <4 x i32> %tmp52, %conv62          ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp63, <4 x i32>* %nominalScore
  %tmp64 = load i8* %encoding_numCodes.addr       ; <i8> [#uses=1]
  %tmp65 = load <2 x i32>** %queryWindowEnd       ; <<2 x i32>*> [#uses=1]
  %conv66 = zext i8 %tmp64 to i32                 ; <i32> [#uses=1]
  %tmp67 = getelementptr <2 x i32>* %tmp65, i32 %conv66 ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp67, <2 x i32>** %queryWindowEnd
  %tmp68 = load i8** %subjectWindowEnd            ; <i8*> [#uses=1]
  %tmp69 = getelementptr i8* %tmp68, i32 1        ; <i8*> [#uses=1]
  store i8* %tmp69, i8** %subjectWindowEnd
  %tmp70 = load <4 x i32>* %count                 ; <<4 x i32>> [#uses=1]
  %tmp71 = add <4 x i32> %tmp70, <i32 1, i32 1, i32 1, i32 1> ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp71, <4 x i32>* %count
  br label %while.cond

while.cond81:                                     ; preds = %if.end140, %while.exit
  %tmp82 = load <2 x i32>** %queryWindowEnd       ; <<2 x i32>*> [#uses=1]
  %structele83 = getelementptr inbounds %struct.PSSMatrixFP* %PSSMatrixFP, i32 0, i32 4 ; <<2 x i32>**> [#uses=1]
  %tmp84 = load <2 x i32>** %structele83          ; <<2 x i32>*> [#uses=1]
  %tmp85 = load %struct.ungappedExtension** %ungappedExtension.addr ; <%struct.ungappedExtension*> [#uses=1]
  %structele86 = getelementptr inbounds %struct.ungappedExtension* %tmp85, i32 0, i32 1 ; <%struct.coordinate*> [#uses=1]
  %structele87 = getelementptr inbounds %struct.coordinate* %structele86, i32 0, i32 0 ; <<4 x i32>*> [#uses=1]
  %tmp88 = load <4 x i32>* %structele87           ; <<4 x i32>> [#uses=1]
  %tmp89 = extractelement <4 x i32> %tmp88, i32 0 ; <i32> [#uses=1]
  %tmp90 = load i8* %encoding_numCodes.addr       ; <i8> [#uses=1]
  %conv91 = zext i8 %tmp90 to i32                 ; <i32> [#uses=1]
  %tmp92 = mul i32 %tmp89, %conv91                ; <i32> [#uses=1]
  %tmp93 = getelementptr <2 x i32>* %tmp84, i32 %tmp92 ; <<2 x i32>*> [#uses=1]
  %cmp94 = icmp ult <2 x i32>* %tmp82, %tmp93     ; <i1> [#uses=1]
  br i1 %cmp94, label %while.body96, label %while.exit95

while.exit95:                                     ; preds = %while.cond81
  %structele145 = getelementptr inbounds %struct.coordinate* %seed, i32 0, i32 0 ; <<4 x i32>*> [#uses=1]
  %tmp146 = load <2 x i32>** %bestQueryPosition   ; <<2 x i32>*> [#uses=1]
  %structele147 = getelementptr inbounds %struct.PSSMatrixFP* %PSSMatrixFP, i32 0, i32 4 ; <<2 x i32>**> [#uses=1]
  %tmp148 = load <2 x i32>** %structele147        ; <<2 x i32>*> [#uses=1]
  %sub.ptr.lhs.cast = ptrtoint <2 x i32>* %tmp146 to i32 ; <i32> [#uses=1]
  %sub.ptr.rhs.cast = ptrtoint <2 x i32>* %tmp148 to i32 ; <i32> [#uses=1]
  %sub.ptr.sub = sub i32 %sub.ptr.lhs.cast, %sub.ptr.rhs.cast ; <i32> [#uses=1]
  %sub.ptr.div = sdiv exact i32 %sub.ptr.sub, 8   ; <i32> [#uses=1]
  %tmp149 = load i8* %encoding_numCodes.addr      ; <i8> [#uses=1]
  %conv150 = zext i8 %tmp149 to i32               ; <i32> [#uses=1]
  %tmp151 = sdiv i32 %sub.ptr.div, %conv150       ; <i32> [#uses=1]
  %tmp152 = add i32 %tmp151, 5                    ; <i32> [#uses=4]
  %conv153 = insertelement <4 x i32> undef, i32 %tmp152, i32 0 ; <<4 x i32>> [#uses=1]
  %conv154 = insertelement <4 x i32> %conv153, i32 %tmp152, i32 1 ; <<4 x i32>> [#uses=1]
  %conv155 = insertelement <4 x i32> %conv154, i32 %tmp152, i32 2 ; <<4 x i32>> [#uses=1]
  %conv156 = insertelement <4 x i32> %conv155, i32 %tmp152, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv156, <4 x i32>* %structele145
  %structele157 = getelementptr inbounds %struct.coordinate* %seed, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp158 = load i8** %bestSubjectPosition        ; <i8*> [#uses=1]
  %tmp159 = getelementptr i8* %tmp158, i32 5      ; <i8*> [#uses=1]
  %tmp160 = load i8** %subject.addr               ; <i8*> [#uses=1]
  %sub.ptr.lhs.cast161 = ptrtoint i8* %tmp159 to i32 ; <i32> [#uses=1]
  %sub.ptr.rhs.cast162 = ptrtoint i8* %tmp160 to i32 ; <i32> [#uses=1]
  %sub.ptr.sub163 = sub i32 %sub.ptr.lhs.cast161, %sub.ptr.rhs.cast162 ; <i32> [#uses=4]
  %conv164 = insertelement <4 x i32> undef, i32 %sub.ptr.sub163, i32 0 ; <<4 x i32>> [#uses=1]
  %conv165 = insertelement <4 x i32> %conv164, i32 %sub.ptr.sub163, i32 1 ; <<4 x i32>> [#uses=1]
  %conv166 = insertelement <4 x i32> %conv165, i32 %sub.ptr.sub163, i32 2 ; <<4 x i32>> [#uses=1]
  %conv167 = insertelement <4 x i32> %conv166, i32 %sub.ptr.sub163, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv167, <4 x i32>* %structele157
  br label %if.end

while.body96:                                     ; preds = %while.cond81
  %tmp97 = load i8* %encoding_numCodes.addr       ; <i8> [#uses=1]
  %tmp98 = load <2 x i32>** %queryWindowEnd       ; <<2 x i32>*> [#uses=1]
  %conv99 = zext i8 %tmp97 to i32                 ; <i32> [#uses=1]
  %tmp100 = getelementptr <2 x i32>* %tmp98, i32 %conv99 ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp100, <2 x i32>** %queryWindowEnd
  %tmp101 = load i8** %subjectWindowEnd           ; <i8*> [#uses=1]
  %tmp102 = getelementptr i8* %tmp101, i32 1      ; <i8*> [#uses=1]
  store i8* %tmp102, i8** %subjectWindowEnd
  %tmp103 = load <4 x i32>* %nominalScore         ; <<4 x i32>> [#uses=1]
  %tmp104 = load <2 x i32>** %queryWindowEnd      ; <<2 x i32>*> [#uses=1]
  %tmp105 = load i8** %subjectWindowEnd           ; <i8*> [#uses=1]
  %tmp106 = load i8* %tmp105                      ; <i8> [#uses=1]
  %conv107 = zext i8 %tmp106 to i32               ; <i32> [#uses=1]
  %arrayidx108 = getelementptr <2 x i32>* %tmp104, i32 %conv107 ; <<2 x i32>*> [#uses=1]
  %tmp109 = load <2 x i32>* %arrayidx108          ; <<2 x i32>> [#uses=1]
  %tmp110 = extractelement <2 x i32> %tmp109, i32 0 ; <i32> [#uses=4]
  %conv111 = insertelement <4 x i32> undef, i32 %tmp110, i32 0 ; <<4 x i32>> [#uses=1]
  %conv112 = insertelement <4 x i32> %conv111, i32 %tmp110, i32 1 ; <<4 x i32>> [#uses=1]
  %conv113 = insertelement <4 x i32> %conv112, i32 %tmp110, i32 2 ; <<4 x i32>> [#uses=1]
  %conv114 = insertelement <4 x i32> %conv113, i32 %tmp110, i32 3 ; <<4 x i32>> [#uses=1]
  %tmp115 = add <4 x i32> %tmp103, %conv114       ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp115, <4 x i32>* %nominalScore
  %tmp116 = load <4 x i32>* %nominalScore         ; <<4 x i32>> [#uses=1]
  %tmp117 = load <2 x i32>** %queryWindowStart    ; <<2 x i32>*> [#uses=1]
  %tmp118 = load i8** %subjectWindowStart         ; <i8*> [#uses=1]
  %tmp119 = load i8* %tmp118                      ; <i8> [#uses=1]
  %conv120 = zext i8 %tmp119 to i32               ; <i32> [#uses=1]
  %arrayidx121 = getelementptr <2 x i32>* %tmp117, i32 %conv120 ; <<2 x i32>*> [#uses=1]
  %tmp122 = load <2 x i32>* %arrayidx121          ; <<2 x i32>> [#uses=1]
  %tmp123 = extractelement <2 x i32> %tmp122, i32 0 ; <i32> [#uses=4]
  %conv124 = insertelement <4 x i32> undef, i32 %tmp123, i32 0 ; <<4 x i32>> [#uses=1]
  %conv125 = insertelement <4 x i32> %conv124, i32 %tmp123, i32 1 ; <<4 x i32>> [#uses=1]
  %conv126 = insertelement <4 x i32> %conv125, i32 %tmp123, i32 2 ; <<4 x i32>> [#uses=1]
  %conv127 = insertelement <4 x i32> %conv126, i32 %tmp123, i32 3 ; <<4 x i32>> [#uses=1]
  %tmp128 = sub <4 x i32> %tmp116, %conv127       ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp128, <4 x i32>* %nominalScore
  %tmp129 = load i8* %encoding_numCodes.addr      ; <i8> [#uses=1]
  %tmp130 = load <2 x i32>** %queryWindowStart    ; <<2 x i32>*> [#uses=1]
  %conv131 = zext i8 %tmp129 to i32               ; <i32> [#uses=1]
  %tmp132 = getelementptr <2 x i32>* %tmp130, i32 %conv131 ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp132, <2 x i32>** %queryWindowStart
  %tmp133 = load i8** %subjectWindowStart         ; <i8*> [#uses=1]
  %tmp134 = getelementptr i8* %tmp133, i32 1      ; <i8*> [#uses=1]
  store i8* %tmp134, i8** %subjectWindowStart
  %tmp135 = load <4 x i32>* %nominalScore         ; <<4 x i32>> [#uses=1]
  %tmp136 = extractelement <4 x i32> %tmp135, i32 0 ; <i32> [#uses=1]
  %tmp137 = load <4 x i32>* %bestSegmentScore     ; <<4 x i32>> [#uses=1]
  %tmp138 = extractelement <4 x i32> %tmp137, i32 0 ; <i32> [#uses=1]
  %cmp139 = icmp sgt i32 %tmp136, %tmp138         ; <i1> [#uses=1]
  br i1 %cmp139, label %if.then141, label %if.end140

if.end140:                                        ; preds = %if.then141, %while.body96
  br label %while.cond81

if.then141:                                       ; preds = %while.body96
  %tmp142 = load <4 x i32>* %nominalScore         ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp142, <4 x i32>* %bestSegmentScore
  %tmp143 = load <2 x i32>** %queryWindowStart    ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp143, <2 x i32>** %bestQueryPosition
  %tmp144 = load i8** %subjectWindowStart         ; <i8*> [#uses=1]
  store i8* %tmp144, i8** %bestSubjectPosition
  br label %if.end140
}

define %struct.ungappedExtension addrspace(1)* @ungappedExtension_twoHitExtendD(i8* %sequenceStart, <4 x i32> %queryOffset, i8* %subjectHit, <4 x i32> %lastHitFP, %struct.PSSMatrixFP* byval %PSSMatrixFP, i8* %subject, <4 x i32>* %sequenceHitEnd, i8 zeroext %encoding_numCodes, <4 x i32> %statistics_ungappedNominalDropoff, <4 x i32> %blast_ungappedNominalTrigger, <4 x i32> %ungappedExtensionsPerThread, %struct.ungappedExtension* %ungappedExtension_extensions, %struct.ungappedExtension* %ungappedExtension_additonal, <4 x i32>* %numOfTriggerExtensions, <4 x i32> %sequenceCount, <4 x i32> %tid) nounwind {
entry:
  %retval = alloca %struct.ungappedExtension addrspace(1)*, align 4 ; <%struct.ungappedExtension addrspace(1)**> [#uses=4]
  %sequenceStart.addr = alloca i8*, align 4       ; <i8**> [#uses=4]
  %queryOffset.addr = alloca <4 x i32>, align 16  ; <<4 x i32>*> [#uses=4]
  %subjectHit.addr = alloca i8*, align 4          ; <i8**> [#uses=6]
  %lastHitFP.addr = alloca <4 x i32>, align 16    ; <<4 x i32>*> [#uses=2]
  %subject.addr = alloca i8*, align 4             ; <i8**> [#uses=5]
  %sequenceHitEnd.addr = alloca <4 x i32>*, align 4 ; <<4 x i32>**> [#uses=3]
  %encoding_numCodes.addr = alloca i8, align 1    ; <i8*> [#uses=8]
  %statistics_ungappedNominalDropoff.addr = alloca <4 x i32>, align 16 ; <<4 x i32>*> [#uses=2]
  %blast_ungappedNominalTrigger.addr = alloca <4 x i32>, align 16 ; <<4 x i32>*> [#uses=2]
  %ungappedExtensionsPerThread.addr = alloca <4 x i32>, align 16 ; <<4 x i32>*> [#uses=2]
  %ungappedExtension_extensions.addr = alloca %struct.ungappedExtension*, align 4 ; <%struct.ungappedExtension**> [#uses=2]
  %ungappedExtension_additonal.addr = alloca %struct.ungappedExtension*, align 4 ; <%struct.ungappedExtension**> [#uses=2]
  %numOfTriggerExtensions.addr = alloca <4 x i32>*, align 4 ; <<4 x i32>**> [#uses=4]
  %sequenceCount.addr = alloca <4 x i32>, align 16 ; <<4 x i32>*> [#uses=2]
  %tid.addr = alloca <4 x i32>, align 16          ; <<4 x i32>*> [#uses=5]
  %queryPosition = alloca <2 x i32>*, align 4     ; <<2 x i32>**> [#uses=14]
  %subjectPosition = alloca i8*, align 4          ; <i8**> [#uses=16]
  %subjectStart = alloca i8*, align 4             ; <i8**> [#uses=6]
  %subjectEnd = alloca i8*, align 4               ; <i8**> [#uses=6]
  %dropoff = alloca <4 x i32>, align 16           ; <<4 x i32>*> [#uses=6]
  %originalDropoff = alloca <4 x i32>, align 16   ; <<4 x i32>*> [#uses=4]
  %ungappedExtension_bestScore = alloca <4 x i32>, align 16 ; <<4 x i32>*> [#uses=10]
  %changeSinceBest = alloca <4 x i32>, align 16   ; <<4 x i32>*> [#uses=16]
  %tmp82 = alloca <4 x i32>, align 16             ; <<4 x i32>*> [#uses=2]
  %tmp162 = alloca <4 x i32>, align 1             ; <<4 x i32>*> [#uses=2]
  %tmp196 = alloca <4 x i32>, align 16            ; <<4 x i32>*> [#uses=9]
  %diagonal = alloca <4 x i32>, align 16          ; <<4 x i32>*> [#uses=3]
  %newUngappedExtension = alloca %struct.ungappedExtension*, align 4 ; <%struct.ungappedExtension**> [#uses=17]
  %tmp268 = alloca %struct.coordinate, align 16   ; <%struct.coordinate*> [#uses=2]
  %tmp271 = alloca %struct.PSSMatrixFP, align 16  ; <%struct.PSSMatrixFP*> [#uses=2]
  store i8* %sequenceStart, i8** %sequenceStart.addr
  store <4 x i32> %queryOffset, <4 x i32>* %queryOffset.addr
  store i8* %subjectHit, i8** %subjectHit.addr
  store <4 x i32> %lastHitFP, <4 x i32>* %lastHitFP.addr
  store i8* %subject, i8** %subject.addr
  store <4 x i32>* %sequenceHitEnd, <4 x i32>** %sequenceHitEnd.addr
  store i8 %encoding_numCodes, i8* %encoding_numCodes.addr
  store <4 x i32> %statistics_ungappedNominalDropoff, <4 x i32>* %statistics_ungappedNominalDropoff.addr
  store <4 x i32> %blast_ungappedNominalTrigger, <4 x i32>* %blast_ungappedNominalTrigger.addr
  store <4 x i32> %ungappedExtensionsPerThread, <4 x i32>* %ungappedExtensionsPerThread.addr
  store %struct.ungappedExtension* %ungappedExtension_extensions, %struct.ungappedExtension** %ungappedExtension_extensions.addr
  store %struct.ungappedExtension* %ungappedExtension_additonal, %struct.ungappedExtension** %ungappedExtension_additonal.addr
  store <4 x i32>* %numOfTriggerExtensions, <4 x i32>** %numOfTriggerExtensions.addr
  store <4 x i32> %sequenceCount, <4 x i32>* %sequenceCount.addr
  store <4 x i32> %tid, <4 x i32>* %tid.addr
  store <4 x i32> zeroinitializer, <4 x i32>* %changeSinceBest
  %tmp = load <4 x i32>* %statistics_ungappedNominalDropoff.addr ; <<4 x i32>> [#uses=1]
  %tmp1 = sub <4 x i32> zeroinitializer, %tmp     ; <<4 x i32>> [#uses=2]
  store <4 x i32> %tmp1, <4 x i32>* %dropoff
  store <4 x i32> %tmp1, <4 x i32>* %originalDropoff
  store <4 x i32> zeroinitializer, <4 x i32>* %ungappedExtension_bestScore
  %structele = getelementptr inbounds %struct.PSSMatrixFP* %PSSMatrixFP, i32 0, i32 4 ; <<2 x i32>**> [#uses=1]
  %tmp2 = load <2 x i32>** %structele             ; <<2 x i32>*> [#uses=1]
  %tmp3 = load <4 x i32>* %queryOffset.addr       ; <<4 x i32>> [#uses=1]
  %tmp4 = extractelement <4 x i32> %tmp3, i32 0   ; <i32> [#uses=1]
  %tmp5 = load i8* %encoding_numCodes.addr        ; <i8> [#uses=1]
  %conv = zext i8 %tmp5 to i32                    ; <i32> [#uses=1]
  %tmp6 = mul i32 %tmp4, %conv                    ; <i32> [#uses=1]
  %tmp7 = getelementptr <2 x i32>* %tmp2, i32 %tmp6 ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp7, <2 x i32>** %queryPosition
  %tmp8 = load i8** %subjectHit.addr              ; <i8*> [#uses=2]
  store i8* %tmp8, i8** %subjectStart
  store i8* %tmp8, i8** %subjectPosition
  br label %while.cond

return:                                           ; preds = %if.else, %if.end203, %if.then67
  %tmp310 = load %struct.ungappedExtension addrspace(1)** %retval ; <%struct.ungappedExtension addrspace(1)*> [#uses=1]
  ret %struct.ungappedExtension addrspace(1)* %tmp310

while.cond:                                       ; preds = %if.end, %entry
  %tmp9 = load <4 x i32>* %changeSinceBest        ; <<4 x i32>> [#uses=1]
  %tmp10 = extractelement <4 x i32> %tmp9, i32 0  ; <i32> [#uses=1]
  %tmp11 = load <4 x i32>* %dropoff               ; <<4 x i32>> [#uses=1]
  %tmp12 = extractelement <4 x i32> %tmp11, i32 0 ; <i32> [#uses=1]
  %cmp = icmp sgt i32 %tmp10, %tmp12              ; <i1> [#uses=1]
  br i1 %cmp, label %while.body, label %while.exit

while.exit:                                       ; preds = %while.cond
  %tmp59 = load i8** %subjectStart                ; <i8*> [#uses=1]
  %tmp60 = getelementptr i8* %tmp59, i32 1        ; <i8*> [#uses=1]
  store i8* %tmp60, i8** %subjectStart
  %tmp61 = load i8** %subjectStart                ; <i8*> [#uses=1]
  %tmp62 = load i8** %sequenceStart.addr          ; <i8*> [#uses=1]
  %sub.ptr.lhs.cast = ptrtoint i8* %tmp61 to i32  ; <i32> [#uses=1]
  %sub.ptr.rhs.cast = ptrtoint i8* %tmp62 to i32  ; <i32> [#uses=1]
  %sub.ptr.sub = sub i32 %sub.ptr.lhs.cast, %sub.ptr.rhs.cast ; <i32> [#uses=1]
  %tmp63 = load <4 x i32>* %lastHitFP.addr        ; <<4 x i32>> [#uses=1]
  %tmp64 = extractelement <4 x i32> %tmp63, i32 0 ; <i32> [#uses=1]
  %cmp65 = icmp ugt i32 %sub.ptr.sub, %tmp64      ; <i1> [#uses=1]
  br i1 %cmp65, label %if.then67, label %if.end66

while.body:                                       ; preds = %while.cond
  %tmp13 = load <4 x i32>* %changeSinceBest       ; <<4 x i32>> [#uses=1]
  %tmp14 = load <2 x i32>** %queryPosition        ; <<2 x i32>*> [#uses=1]
  %tmp15 = load i8** %subjectPosition             ; <i8*> [#uses=1]
  %tmp16 = load i8* %tmp15                        ; <i8> [#uses=1]
  %conv17 = zext i8 %tmp16 to i32                 ; <i32> [#uses=1]
  %arrayidx = getelementptr <2 x i32>* %tmp14, i32 %conv17 ; <<2 x i32>*> [#uses=1]
  %tmp18 = load <2 x i32>* %arrayidx              ; <<2 x i32>> [#uses=1]
  %tmp19 = extractelement <2 x i32> %tmp18, i32 0 ; <i32> [#uses=4]
  %conv20 = insertelement <4 x i32> undef, i32 %tmp19, i32 0 ; <<4 x i32>> [#uses=1]
  %conv21 = insertelement <4 x i32> %conv20, i32 %tmp19, i32 1 ; <<4 x i32>> [#uses=1]
  %conv22 = insertelement <4 x i32> %conv21, i32 %tmp19, i32 2 ; <<4 x i32>> [#uses=1]
  %conv23 = insertelement <4 x i32> %conv22, i32 %tmp19, i32 3 ; <<4 x i32>> [#uses=1]
  %tmp24 = add <4 x i32> %tmp13, %conv23          ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp24, <4 x i32>* %changeSinceBest
  %tmp25 = load <4 x i32>* %changeSinceBest       ; <<4 x i32>> [#uses=1]
  %tmp26 = extractelement <4 x i32> %tmp25, i32 0 ; <i32> [#uses=1]
  %cmp27 = icmp sgt i32 %tmp26, 0                 ; <i1> [#uses=1]
  br i1 %cmp27, label %if.then, label %if.end

if.end:                                           ; preds = %do.exit, %while.body
  %tmp52 = load <2 x i32>** %queryPosition        ; <<2 x i32>*> [#uses=1]
  %tmp53 = load i8* %encoding_numCodes.addr       ; <i8> [#uses=1]
  %conv54 = zext i8 %tmp53 to i32                 ; <i32> [#uses=1]
  %sub.ptr.neg55 = sub i32 0, %conv54             ; <i32> [#uses=1]
  %tmp56 = getelementptr <2 x i32>* %tmp52, i32 %sub.ptr.neg55 ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp56, <2 x i32>** %queryPosition
  %tmp57 = load i8** %subjectPosition             ; <i8*> [#uses=1]
  %tmp58 = getelementptr i8* %tmp57, i32 -1       ; <i8*> [#uses=1]
  store i8* %tmp58, i8** %subjectPosition
  br label %while.cond

if.then:                                          ; preds = %while.body
  br label %do.body

do.exit:                                          ; preds = %do.cond
  %tmp51 = load i8** %subjectPosition             ; <i8*> [#uses=1]
  store i8* %tmp51, i8** %subjectStart
  br label %if.end

do.body:                                          ; preds = %do.cond, %if.then
  %tmp28 = load <4 x i32>* %ungappedExtension_bestScore ; <<4 x i32>> [#uses=1]
  %tmp29 = load <4 x i32>* %changeSinceBest       ; <<4 x i32>> [#uses=1]
  %tmp30 = add <4 x i32> %tmp28, %tmp29           ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp30, <4 x i32>* %ungappedExtension_bestScore
  %tmp31 = load <2 x i32>** %queryPosition        ; <<2 x i32>*> [#uses=1]
  %tmp32 = load i8* %encoding_numCodes.addr       ; <i8> [#uses=1]
  %conv33 = zext i8 %tmp32 to i32                 ; <i32> [#uses=1]
  %sub.ptr.neg = sub i32 0, %conv33               ; <i32> [#uses=1]
  %tmp34 = getelementptr <2 x i32>* %tmp31, i32 %sub.ptr.neg ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp34, <2 x i32>** %queryPosition
  %tmp35 = load i8** %subjectPosition             ; <i8*> [#uses=1]
  %tmp36 = getelementptr i8* %tmp35, i32 -1       ; <i8*> [#uses=1]
  store i8* %tmp36, i8** %subjectPosition
  %tmp37 = load <2 x i32>** %queryPosition        ; <<2 x i32>*> [#uses=1]
  %tmp38 = load i8** %subjectPosition             ; <i8*> [#uses=1]
  %tmp39 = load i8* %tmp38                        ; <i8> [#uses=1]
  %conv40 = zext i8 %tmp39 to i32                 ; <i32> [#uses=1]
  %arrayidx41 = getelementptr <2 x i32>* %tmp37, i32 %conv40 ; <<2 x i32>*> [#uses=1]
  %tmp42 = load <2 x i32>* %arrayidx41            ; <<2 x i32>> [#uses=1]
  %tmp43 = extractelement <2 x i32> %tmp42, i32 0 ; <i32> [#uses=4]
  %conv44 = insertelement <4 x i32> undef, i32 %tmp43, i32 0 ; <<4 x i32>> [#uses=1]
  %conv45 = insertelement <4 x i32> %conv44, i32 %tmp43, i32 1 ; <<4 x i32>> [#uses=1]
  %conv46 = insertelement <4 x i32> %conv45, i32 %tmp43, i32 2 ; <<4 x i32>> [#uses=1]
  %conv47 = insertelement <4 x i32> %conv46, i32 %tmp43, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv47, <4 x i32>* %changeSinceBest
  br label %do.cond

do.cond:                                          ; preds = %do.body
  %tmp48 = load <4 x i32>* %changeSinceBest       ; <<4 x i32>> [#uses=1]
  %tmp49 = extractelement <4 x i32> %tmp48, i32 0 ; <i32> [#uses=1]
  %cmp50 = icmp sgt i32 %tmp49, 0                 ; <i1> [#uses=1]
  br i1 %cmp50, label %do.body, label %do.exit

if.end66:                                         ; preds = %while.exit
  %structele78 = getelementptr inbounds %struct.PSSMatrixFP* %PSSMatrixFP, i32 0, i32 4 ; <<2 x i32>**> [#uses=1]
  %tmp79 = load <2 x i32>** %structele78          ; <<2 x i32>*> [#uses=1]
  %tmp80 = load <4 x i32>* %queryOffset.addr      ; <<4 x i32>> [#uses=1]
  %tmp81 = add <4 x i32> %tmp80, <i32 1, i32 1, i32 1, i32 1> ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp81, <4 x i32>* %tmp82
  %tmp83 = load <4 x i32>* %tmp82                 ; <<4 x i32>> [#uses=1]
  %tmp84 = extractelement <4 x i32> %tmp83, i32 0 ; <i32> [#uses=1]
  %tmp85 = load i8* %encoding_numCodes.addr       ; <i8> [#uses=1]
  %conv86 = zext i8 %tmp85 to i32                 ; <i32> [#uses=1]
  %tmp87 = mul i32 %tmp84, %conv86                ; <i32> [#uses=1]
  %tmp88 = getelementptr <2 x i32>* %tmp79, i32 %tmp87 ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp88, <2 x i32>** %queryPosition
  %tmp89 = load i8** %subjectHit.addr             ; <i8*> [#uses=1]
  %tmp90 = getelementptr i8* %tmp89, i32 1        ; <i8*> [#uses=1]
  store i8* %tmp90, i8** %subjectPosition
  %tmp91 = load i8** %subjectHit.addr             ; <i8*> [#uses=1]
  store i8* %tmp91, i8** %subjectEnd
  store <4 x i32> zeroinitializer, <4 x i32>* %changeSinceBest
  %tmp92 = load <4 x i32>* %ungappedExtension_bestScore ; <<4 x i32>> [#uses=1]
  %tmp93 = extractelement <4 x i32> %tmp92, i32 0 ; <i32> [#uses=1]
  %tmp94 = sub i32 0, %tmp93                      ; <i32> [#uses=1]
  %tmp95 = load <4 x i32>* %originalDropoff       ; <<4 x i32>> [#uses=1]
  %tmp96 = extractelement <4 x i32> %tmp95, i32 0 ; <i32> [#uses=1]
  %cmp97 = icmp sgt i32 %tmp94, %tmp96            ; <i1> [#uses=1]
  br i1 %cmp97, label %if.then99, label %if.end98

if.then67:                                        ; preds = %while.exit
  %tmp68 = load <4 x i32>** %sequenceHitEnd.addr  ; <<4 x i32>*> [#uses=1]
  %tmp69 = load i8** %subjectHit.addr             ; <i8*> [#uses=1]
  %tmp70 = load i8** %sequenceStart.addr          ; <i8*> [#uses=1]
  %sub.ptr.lhs.cast71 = ptrtoint i8* %tmp69 to i32 ; <i32> [#uses=1]
  %sub.ptr.rhs.cast72 = ptrtoint i8* %tmp70 to i32 ; <i32> [#uses=1]
  %sub.ptr.sub73 = sub i32 %sub.ptr.lhs.cast71, %sub.ptr.rhs.cast72 ; <i32> [#uses=4]
  %conv74 = insertelement <4 x i32> undef, i32 %sub.ptr.sub73, i32 0 ; <<4 x i32>> [#uses=1]
  %conv75 = insertelement <4 x i32> %conv74, i32 %sub.ptr.sub73, i32 1 ; <<4 x i32>> [#uses=1]
  %conv76 = insertelement <4 x i32> %conv75, i32 %sub.ptr.sub73, i32 2 ; <<4 x i32>> [#uses=1]
  %conv77 = insertelement <4 x i32> %conv76, i32 %sub.ptr.sub73, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv77, <4 x i32>* %tmp68
  store %struct.ungappedExtension addrspace(1)* null, %struct.ungappedExtension addrspace(1)** %retval
  br label %return

if.end98:                                         ; preds = %if.then99, %if.end66
  br label %while.cond102

if.then99:                                        ; preds = %if.end66
  %tmp100 = load <4 x i32>* %ungappedExtension_bestScore ; <<4 x i32>> [#uses=1]
  %tmp101 = sub <4 x i32> zeroinitializer, %tmp100 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp101, <4 x i32>* %dropoff
  br label %if.end98

while.cond102:                                    ; preds = %if.end126, %if.end98
  %tmp103 = load <4 x i32>* %changeSinceBest      ; <<4 x i32>> [#uses=1]
  %tmp104 = extractelement <4 x i32> %tmp103, i32 0 ; <i32> [#uses=1]
  %tmp105 = load <4 x i32>* %dropoff              ; <<4 x i32>> [#uses=1]
  %tmp106 = extractelement <4 x i32> %tmp105, i32 0 ; <i32> [#uses=1]
  %cmp107 = icmp sgt i32 %tmp104, %tmp106         ; <i1> [#uses=1]
  br i1 %cmp107, label %while.body109, label %while.exit108

while.exit108:                                    ; preds = %while.cond102
  %tmp177 = load i8** %subjectEnd                 ; <i8*> [#uses=1]
  %tmp178 = getelementptr i8* %tmp177, i32 -1     ; <i8*> [#uses=1]
  store i8* %tmp178, i8** %subjectEnd
  %tmp179 = load <4 x i32>** %sequenceHitEnd.addr ; <<4 x i32>*> [#uses=1]
  %tmp180 = load i8** %subjectEnd                 ; <i8*> [#uses=1]
  %tmp181 = load i8** %sequenceStart.addr         ; <i8*> [#uses=1]
  %sub.ptr.lhs.cast182 = ptrtoint i8* %tmp180 to i32 ; <i32> [#uses=1]
  %sub.ptr.rhs.cast183 = ptrtoint i8* %tmp181 to i32 ; <i32> [#uses=1]
  %sub.ptr.sub184 = sub i32 %sub.ptr.lhs.cast182, %sub.ptr.rhs.cast183 ; <i32> [#uses=4]
  %conv185 = insertelement <4 x i32> undef, i32 %sub.ptr.sub184, i32 0 ; <<4 x i32>> [#uses=1]
  %conv186 = insertelement <4 x i32> %conv185, i32 %sub.ptr.sub184, i32 1 ; <<4 x i32>> [#uses=1]
  %conv187 = insertelement <4 x i32> %conv186, i32 %sub.ptr.sub184, i32 2 ; <<4 x i32>> [#uses=1]
  %conv188 = insertelement <4 x i32> %conv187, i32 %sub.ptr.sub184, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv188, <4 x i32>* %tmp179
  %tmp189 = load <4 x i32>* %ungappedExtension_bestScore ; <<4 x i32>> [#uses=1]
  %tmp190 = extractelement <4 x i32> %tmp189, i32 0 ; <i32> [#uses=1]
  %tmp191 = load <4 x i32>* %blast_ungappedNominalTrigger.addr ; <<4 x i32>> [#uses=1]
  %tmp192 = extractelement <4 x i32> %tmp191, i32 0 ; <i32> [#uses=1]
  %cmp193 = icmp sge i32 %tmp190, %tmp192         ; <i1> [#uses=1]
  br i1 %cmp193, label %if.then195, label %if.else

while.body109:                                    ; preds = %while.cond102
  %tmp110 = load <4 x i32>* %changeSinceBest      ; <<4 x i32>> [#uses=1]
  %tmp111 = load <2 x i32>** %queryPosition       ; <<2 x i32>*> [#uses=1]
  %tmp112 = load i8** %subjectPosition            ; <i8*> [#uses=1]
  %tmp113 = load i8* %tmp112                      ; <i8> [#uses=1]
  %conv114 = zext i8 %tmp113 to i32               ; <i32> [#uses=1]
  %arrayidx115 = getelementptr <2 x i32>* %tmp111, i32 %conv114 ; <<2 x i32>*> [#uses=1]
  %tmp116 = load <2 x i32>* %arrayidx115          ; <<2 x i32>> [#uses=1]
  %tmp117 = extractelement <2 x i32> %tmp116, i32 0 ; <i32> [#uses=4]
  %conv118 = insertelement <4 x i32> undef, i32 %tmp117, i32 0 ; <<4 x i32>> [#uses=1]
  %conv119 = insertelement <4 x i32> %conv118, i32 %tmp117, i32 1 ; <<4 x i32>> [#uses=1]
  %conv120 = insertelement <4 x i32> %conv119, i32 %tmp117, i32 2 ; <<4 x i32>> [#uses=1]
  %conv121 = insertelement <4 x i32> %conv120, i32 %tmp117, i32 3 ; <<4 x i32>> [#uses=1]
  %tmp122 = add <4 x i32> %tmp110, %conv121       ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp122, <4 x i32>* %changeSinceBest
  %tmp123 = load <4 x i32>* %changeSinceBest      ; <<4 x i32>> [#uses=1]
  %tmp124 = extractelement <4 x i32> %tmp123, i32 0 ; <i32> [#uses=1]
  %cmp125 = icmp sgt i32 %tmp124, 0               ; <i1> [#uses=1]
  br i1 %cmp125, label %if.then127, label %if.end126

if.end126:                                        ; preds = %if.end168, %while.body109
  %tmp171 = load <2 x i32>** %queryPosition       ; <<2 x i32>*> [#uses=1]
  %tmp172 = load i8* %encoding_numCodes.addr      ; <i8> [#uses=1]
  %conv173 = zext i8 %tmp172 to i32               ; <i32> [#uses=1]
  %tmp174 = getelementptr <2 x i32>* %tmp171, i32 %conv173 ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp174, <2 x i32>** %queryPosition
  %tmp175 = load i8** %subjectPosition            ; <i8*> [#uses=1]
  %tmp176 = getelementptr i8* %tmp175, i32 1      ; <i8*> [#uses=1]
  store i8* %tmp176, i8** %subjectPosition
  br label %while.cond102

if.then127:                                       ; preds = %while.body109
  br label %do.body129

do.exit128:                                       ; preds = %do.cond130
  %tmp154 = load i8** %subjectPosition            ; <i8*> [#uses=1]
  store i8* %tmp154, i8** %subjectEnd
  %tmp155 = load <4 x i32>* %ungappedExtension_bestScore ; <<4 x i32>> [#uses=1]
  %tmp156 = extractelement <4 x i32> %tmp155, i32 0 ; <i32> [#uses=1]
  %tmp157 = sub i32 0, %tmp156                    ; <i32> [#uses=4]
  %conv158 = insertelement <4 x i32> undef, i32 %tmp157, i32 0 ; <<4 x i32>> [#uses=1]
  %conv159 = insertelement <4 x i32> %conv158, i32 %tmp157, i32 1 ; <<4 x i32>> [#uses=1]
  %conv160 = insertelement <4 x i32> %conv159, i32 %tmp157, i32 2 ; <<4 x i32>> [#uses=1]
  %conv161 = insertelement <4 x i32> %conv160, i32 %tmp157, i32 3 ; <<4 x i32>> [#uses=2]
  store <4 x i32> %conv161, <4 x i32>* %dropoff
  store <4 x i32> %conv161, <4 x i32>* %tmp162
  %tmp163 = load <4 x i32>* %tmp162               ; <<4 x i32>> [#uses=1]
  %tmp164 = extractelement <4 x i32> %tmp163, i32 0 ; <i32> [#uses=1]
  %tmp165 = load <4 x i32>* %originalDropoff      ; <<4 x i32>> [#uses=1]
  %tmp166 = extractelement <4 x i32> %tmp165, i32 0 ; <i32> [#uses=1]
  %cmp167 = icmp slt i32 %tmp164, %tmp166         ; <i1> [#uses=1]
  br i1 %cmp167, label %if.then169, label %if.end168

do.body129:                                       ; preds = %do.cond130, %if.then127
  %tmp131 = load <4 x i32>* %ungappedExtension_bestScore ; <<4 x i32>> [#uses=1]
  %tmp132 = load <4 x i32>* %changeSinceBest      ; <<4 x i32>> [#uses=1]
  %tmp133 = add <4 x i32> %tmp131, %tmp132        ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp133, <4 x i32>* %ungappedExtension_bestScore
  %tmp134 = load <2 x i32>** %queryPosition       ; <<2 x i32>*> [#uses=1]
  %tmp135 = load i8* %encoding_numCodes.addr      ; <i8> [#uses=1]
  %conv136 = zext i8 %tmp135 to i32               ; <i32> [#uses=1]
  %tmp137 = getelementptr <2 x i32>* %tmp134, i32 %conv136 ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp137, <2 x i32>** %queryPosition
  %tmp138 = load i8** %subjectPosition            ; <i8*> [#uses=1]
  %tmp139 = getelementptr i8* %tmp138, i32 1      ; <i8*> [#uses=1]
  store i8* %tmp139, i8** %subjectPosition
  %tmp140 = load <2 x i32>** %queryPosition       ; <<2 x i32>*> [#uses=1]
  %tmp141 = load i8** %subjectPosition            ; <i8*> [#uses=1]
  %tmp142 = load i8* %tmp141                      ; <i8> [#uses=1]
  %conv143 = zext i8 %tmp142 to i32               ; <i32> [#uses=1]
  %arrayidx144 = getelementptr <2 x i32>* %tmp140, i32 %conv143 ; <<2 x i32>*> [#uses=1]
  %tmp145 = load <2 x i32>* %arrayidx144          ; <<2 x i32>> [#uses=1]
  %tmp146 = extractelement <2 x i32> %tmp145, i32 0 ; <i32> [#uses=4]
  %conv147 = insertelement <4 x i32> undef, i32 %tmp146, i32 0 ; <<4 x i32>> [#uses=1]
  %conv148 = insertelement <4 x i32> %conv147, i32 %tmp146, i32 1 ; <<4 x i32>> [#uses=1]
  %conv149 = insertelement <4 x i32> %conv148, i32 %tmp146, i32 2 ; <<4 x i32>> [#uses=1]
  %conv150 = insertelement <4 x i32> %conv149, i32 %tmp146, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv150, <4 x i32>* %changeSinceBest
  br label %do.cond130

do.cond130:                                       ; preds = %do.body129
  %tmp151 = load <4 x i32>* %changeSinceBest      ; <<4 x i32>> [#uses=1]
  %tmp152 = extractelement <4 x i32> %tmp151, i32 0 ; <i32> [#uses=1]
  %cmp153 = icmp sgt i32 %tmp152, 0               ; <i1> [#uses=1]
  br i1 %cmp153, label %do.body129, label %do.exit128

if.end168:                                        ; preds = %if.then169, %do.exit128
  br label %if.end126

if.then169:                                       ; preds = %do.exit128
  %tmp170 = load <4 x i32>* %originalDropoff      ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp170, <4 x i32>* %dropoff
  br label %if.end168

if.then195:                                       ; preds = %while.exit108
  store %struct.ungappedExtension* null, %struct.ungappedExtension** %newUngappedExtension
  %tmp197 = load <4 x i32>** %numOfTriggerExtensions.addr ; <<4 x i32>*> [#uses=1]
  %tmp198 = load <4 x i32>* %tmp197               ; <<4 x i32>> [#uses=1]
  %tmp199 = extractelement <4 x i32> %tmp198, i32 0 ; <i32> [#uses=1]
  %tmp200 = load <4 x i32>* %ungappedExtensionsPerThread.addr ; <<4 x i32>> [#uses=1]
  %tmp201 = extractelement <4 x i32> %tmp200, i32 0 ; <i32> [#uses=1]
  %cmp202 = icmp uge i32 %tmp199, %tmp201         ; <i1> [#uses=1]
  br i1 %cmp202, label %if.then204, label %if.else205

if.else:                                          ; preds = %while.exit108
  store %struct.ungappedExtension addrspace(1)* null, %struct.ungappedExtension addrspace(1)** %retval
  br label %return

if.end203:                                        ; preds = %if.else205, %if.then204
  %tmp213 = load i8** %subjectHit.addr            ; <i8*> [#uses=1]
  %tmp214 = load i8** %subject.addr               ; <i8*> [#uses=1]
  %sub.ptr.lhs.cast215 = ptrtoint i8* %tmp213 to i32 ; <i32> [#uses=1]
  %sub.ptr.rhs.cast216 = ptrtoint i8* %tmp214 to i32 ; <i32> [#uses=1]
  %sub.ptr.sub217 = sub i32 %sub.ptr.lhs.cast215, %sub.ptr.rhs.cast216 ; <i32> [#uses=4]
  %conv218 = insertelement <4 x i32> undef, i32 %sub.ptr.sub217, i32 0 ; <<4 x i32>> [#uses=1]
  %conv219 = insertelement <4 x i32> %conv218, i32 %sub.ptr.sub217, i32 1 ; <<4 x i32>> [#uses=1]
  %conv220 = insertelement <4 x i32> %conv219, i32 %sub.ptr.sub217, i32 2 ; <<4 x i32>> [#uses=1]
  %conv221 = insertelement <4 x i32> %conv220, i32 %sub.ptr.sub217, i32 3 ; <<4 x i32>> [#uses=1]
  %tmp222 = load <4 x i32>* %queryOffset.addr     ; <<4 x i32>> [#uses=1]
  %tmp223 = sub <4 x i32> %conv221, %tmp222       ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp223, <4 x i32>* %diagonal
  %tmp224 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele225 = getelementptr inbounds %struct.ungappedExtension* %tmp224, i32 0, i32 0 ; <%struct.coordinate*> [#uses=1]
  %structele226 = getelementptr inbounds %struct.coordinate* %structele225, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp227 = load i8** %subjectStart               ; <i8*> [#uses=1]
  %tmp228 = load i8** %subject.addr               ; <i8*> [#uses=1]
  %sub.ptr.lhs.cast229 = ptrtoint i8* %tmp227 to i32 ; <i32> [#uses=1]
  %sub.ptr.rhs.cast230 = ptrtoint i8* %tmp228 to i32 ; <i32> [#uses=1]
  %sub.ptr.sub231 = sub i32 %sub.ptr.lhs.cast229, %sub.ptr.rhs.cast230 ; <i32> [#uses=4]
  %conv232 = insertelement <4 x i32> undef, i32 %sub.ptr.sub231, i32 0 ; <<4 x i32>> [#uses=1]
  %conv233 = insertelement <4 x i32> %conv232, i32 %sub.ptr.sub231, i32 1 ; <<4 x i32>> [#uses=1]
  %conv234 = insertelement <4 x i32> %conv233, i32 %sub.ptr.sub231, i32 2 ; <<4 x i32>> [#uses=1]
  %conv235 = insertelement <4 x i32> %conv234, i32 %sub.ptr.sub231, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv235, <4 x i32>* %structele226
  %tmp236 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele237 = getelementptr inbounds %struct.ungappedExtension* %tmp236, i32 0, i32 1 ; <%struct.coordinate*> [#uses=1]
  %structele238 = getelementptr inbounds %struct.coordinate* %structele237, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp239 = load i8** %subjectEnd                 ; <i8*> [#uses=1]
  %tmp240 = load i8** %subject.addr               ; <i8*> [#uses=1]
  %sub.ptr.lhs.cast241 = ptrtoint i8* %tmp239 to i32 ; <i32> [#uses=1]
  %sub.ptr.rhs.cast242 = ptrtoint i8* %tmp240 to i32 ; <i32> [#uses=1]
  %sub.ptr.sub243 = sub i32 %sub.ptr.lhs.cast241, %sub.ptr.rhs.cast242 ; <i32> [#uses=4]
  %conv244 = insertelement <4 x i32> undef, i32 %sub.ptr.sub243, i32 0 ; <<4 x i32>> [#uses=1]
  %conv245 = insertelement <4 x i32> %conv244, i32 %sub.ptr.sub243, i32 1 ; <<4 x i32>> [#uses=1]
  %conv246 = insertelement <4 x i32> %conv245, i32 %sub.ptr.sub243, i32 2 ; <<4 x i32>> [#uses=1]
  %conv247 = insertelement <4 x i32> %conv246, i32 %sub.ptr.sub243, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv247, <4 x i32>* %structele238
  %tmp248 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele249 = getelementptr inbounds %struct.ungappedExtension* %tmp248, i32 0, i32 0 ; <%struct.coordinate*> [#uses=1]
  %structele250 = getelementptr inbounds %struct.coordinate* %structele249, i32 0, i32 0 ; <<4 x i32>*> [#uses=1]
  %tmp251 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele252 = getelementptr inbounds %struct.ungappedExtension* %tmp251, i32 0, i32 0 ; <%struct.coordinate*> [#uses=1]
  %structele253 = getelementptr inbounds %struct.coordinate* %structele252, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp254 = load <4 x i32>* %structele253         ; <<4 x i32>> [#uses=1]
  %tmp255 = load <4 x i32>* %diagonal             ; <<4 x i32>> [#uses=1]
  %tmp256 = sub <4 x i32> %tmp254, %tmp255        ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp256, <4 x i32>* %structele250
  %tmp257 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele258 = getelementptr inbounds %struct.ungappedExtension* %tmp257, i32 0, i32 1 ; <%struct.coordinate*> [#uses=1]
  %structele259 = getelementptr inbounds %struct.coordinate* %structele258, i32 0, i32 0 ; <<4 x i32>*> [#uses=1]
  %tmp260 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele261 = getelementptr inbounds %struct.ungappedExtension* %tmp260, i32 0, i32 1 ; <%struct.coordinate*> [#uses=1]
  %structele262 = getelementptr inbounds %struct.coordinate* %structele261, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp263 = load <4 x i32>* %structele262         ; <<4 x i32>> [#uses=1]
  %tmp264 = load <4 x i32>* %diagonal             ; <<4 x i32>> [#uses=1]
  %tmp265 = sub <4 x i32> %tmp263, %tmp264        ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp265, <4 x i32>* %structele259
  %tmp266 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele267 = getelementptr inbounds %struct.ungappedExtension* %tmp266, i32 0, i32 2 ; <%struct.coordinate*> [#uses=1]
  %tmp269 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %tmp270 = load %struct.PSSMatrixFP* %PSSMatrixFP ; <%struct.PSSMatrixFP> [#uses=1]
  store %struct.PSSMatrixFP %tmp270, %struct.PSSMatrixFP* %tmp271
  %tmp272 = load i8** %subject.addr               ; <i8*> [#uses=1]
  %tmp273 = load i8* %encoding_numCodes.addr      ; <i8> [#uses=1]
  call void @ungappedExtension_findProteinSeed(%struct.coordinate* noalias sret %tmp268, %struct.ungappedExtension* %tmp269, %struct.PSSMatrixFP* byval %tmp271, i8* %tmp272, i8 zeroext %tmp273) nounwind
  %tmp274 = load %struct.coordinate* %tmp268      ; <%struct.coordinate> [#uses=1]
  store %struct.coordinate %tmp274, %struct.coordinate* %structele267
  %tmp275 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele276 = getelementptr inbounds %struct.ungappedExtension* %tmp275, i32 0, i32 10 ; <%struct.ungappedExtension**> [#uses=1]
  store %struct.ungappedExtension* null, %struct.ungappedExtension** %structele276
  %tmp277 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele278 = getelementptr inbounds %struct.ungappedExtension* %tmp277, i32 0, i32 3 ; <<4 x i32>*> [#uses=1]
  %tmp279 = load <4 x i32>* %ungappedExtension_bestScore ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp279, <4 x i32>* %structele278
  %tmp280 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele281 = getelementptr inbounds %struct.ungappedExtension* %tmp280, i32 0, i32 6 ; <i8*> [#uses=1]
  store i8 1, i8* %structele281
  %tmp282 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele283 = getelementptr inbounds %struct.ungappedExtension* %tmp282, i32 0, i32 4 ; <<4 x i32>*> [#uses=1]
  %tmp284 = load <4 x i32>* %sequenceCount.addr   ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp284, <4 x i32>* %structele283
  %tmp285 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %structele286 = getelementptr inbounds %struct.ungappedExtension* %tmp285, i32 0, i32 5 ; <<4 x i32>*> [#uses=1]
  %tmp287 = load <4 x i32>* %tid.addr             ; <<4 x i32>> [#uses=1]
  %tmp288 = extractelement <4 x i32> %tmp287, i32 0 ; <i32> [#uses=1]
  %tmp289 = load <4 x i32>* %tmp196               ; <<4 x i32>> [#uses=1]
  %tmp290 = insertelement <4 x i32> %tmp289, i32 %tmp288, i32 0 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp290, <4 x i32>* %tmp196
  %tmp291 = load <4 x i32>* %tid.addr             ; <<4 x i32>> [#uses=1]
  %tmp292 = extractelement <4 x i32> %tmp291, i32 1 ; <i32> [#uses=1]
  %tmp293 = load <4 x i32>* %tmp196               ; <<4 x i32>> [#uses=1]
  %tmp294 = insertelement <4 x i32> %tmp293, i32 %tmp292, i32 1 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp294, <4 x i32>* %tmp196
  %tmp295 = load <4 x i32>* %tid.addr             ; <<4 x i32>> [#uses=1]
  %tmp296 = extractelement <4 x i32> %tmp295, i32 2 ; <i32> [#uses=1]
  %tmp297 = load <4 x i32>* %tmp196               ; <<4 x i32>> [#uses=1]
  %tmp298 = insertelement <4 x i32> %tmp297, i32 %tmp296, i32 2 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp298, <4 x i32>* %tmp196
  %tmp299 = load <4 x i32>* %tid.addr             ; <<4 x i32>> [#uses=1]
  %tmp300 = extractelement <4 x i32> %tmp299, i32 3 ; <i32> [#uses=1]
  %tmp301 = load <4 x i32>* %tmp196               ; <<4 x i32>> [#uses=1]
  %tmp302 = insertelement <4 x i32> %tmp301, i32 %tmp300, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp302, <4 x i32>* %tmp196
  %tmp303 = load <4 x i32>* %tmp196               ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp303, <4 x i32>* %structele286
  %tmp304 = load <4 x i32>** %numOfTriggerExtensions.addr ; <<4 x i32>*> [#uses=2]
  %tmp305 = load <4 x i32>* %tmp304               ; <<4 x i32>> [#uses=1]
  %tmp306 = add <4 x i32> %tmp305, <i32 1, i32 1, i32 1, i32 1> ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp306, <4 x i32>* %tmp304
  %tmp307 = load %struct.ungappedExtension** %newUngappedExtension ; <%struct.ungappedExtension*> [#uses=1]
  %conv308 = ptrtoint %struct.ungappedExtension* %tmp307 to i32 ; <i32> [#uses=1]
  %conv309 = inttoptr i32 %conv308 to %struct.ungappedExtension addrspace(1)* ; <%struct.ungappedExtension addrspace(1)*> [#uses=1]
  store %struct.ungappedExtension addrspace(1)* %conv309, %struct.ungappedExtension addrspace(1)** %retval
  br label %return

if.then204:                                       ; preds = %if.then195
  %tmp206 = load %struct.ungappedExtension** %ungappedExtension_additonal.addr ; <%struct.ungappedExtension*> [#uses=1]
  %call = call i32 @__atomic_add_gi32(i32 addrspace(1)* inttoptr (i32 ptrtoint (<4 x i32>* @global_numAdditionalTriggerExtensions to i32) to i32 addrspace(1)*), i32 1) nounwind ; <i32> [#uses=1]
  %tmp207 = getelementptr %struct.ungappedExtension* %tmp206, i32 %call ; <%struct.ungappedExtension*> [#uses=1]
  store %struct.ungappedExtension* %tmp207, %struct.ungappedExtension** %newUngappedExtension
  br label %if.end203

if.else205:                                       ; preds = %if.then195
  %tmp208 = load %struct.ungappedExtension** %ungappedExtension_extensions.addr ; <%struct.ungappedExtension*> [#uses=1]
  %tmp209 = load <4 x i32>** %numOfTriggerExtensions.addr ; <<4 x i32>*> [#uses=1]
  %tmp210 = load <4 x i32>* %tmp209               ; <<4 x i32>> [#uses=1]
  %tmp211 = extractelement <4 x i32> %tmp210, i32 0 ; <i32> [#uses=1]
  %tmp212 = getelementptr %struct.ungappedExtension* %tmp208, i32 %tmp211 ; <%struct.ungappedExtension*> [#uses=1]
  store %struct.ungappedExtension* %tmp212, %struct.ungappedExtension** %newUngappedExtension
  br label %if.end203
}

declare i32 @__atomic_add_gi32(i32 addrspace(1)*, i32) nounwind

define void @__OpenCL_search_protein2hitKernel_kernel(%struct.PSSMatrixFP addrspace(1)* %PSSMatrixFP, <2 x i32> addrspace(1)* %matrixBody, %struct.sequenceDataFP addrspace(1)* %sequenceDataFP, i8 addrspace(1)* %sequence, %struct.parameters addrspace(1)* %parametersFP, %struct.coordinate addrspace(1)* %wordLookupDFA_groupsFP, i8 addrspace(1)* %wordLookupDFAFP, <4 x i32> addrspace(1)* %blast_numUngappedExtensions, <4 x i32> addrspace(1)* %blast_numTriggerExtensions, <4 x i32> addrspace(1)* %blast_numHits, <4 x i32> addrspace(1)* %hitMatrix_furthestp, <4 x i32> addrspace(1)* %hitMatrix_offsetp, %struct.ungappedExtension addrspace(1)* %ungappedExtension_extensionsp, <4 x i32> %nTotalSequenceNum) nounwind {
entry:
  %PSSMatrixFP.addr = alloca %struct.PSSMatrixFP addrspace(1)*, align 4 ; <%struct.PSSMatrixFP addrspace(1)**> [#uses=4]
  %matrixBody.addr = alloca <2 x i32> addrspace(1)*, align 4 ; <<2 x i32> addrspace(1)**> [#uses=2]
  %sequenceDataFP.addr = alloca %struct.sequenceDataFP addrspace(1)*, align 4 ; <%struct.sequenceDataFP addrspace(1)**> [#uses=4]
  %sequence.addr = alloca i8 addrspace(1)*, align 4 ; <i8 addrspace(1)**> [#uses=5]
  %parametersFP.addr = alloca %struct.parameters addrspace(1)*, align 4 ; <%struct.parameters addrspace(1)**> [#uses=15]
  %wordLookupDFA_groupsFP.addr = alloca %struct.coordinate addrspace(1)*, align 4 ; <%struct.coordinate addrspace(1)**> [#uses=6]
  %wordLookupDFAFP.addr = alloca i8 addrspace(1)*, align 4 ; <i8 addrspace(1)**> [#uses=3]
  %blast_numUngappedExtensions.addr = alloca <4 x i32> addrspace(1)*, align 4 ; <<4 x i32> addrspace(1)**> [#uses=2]
  %blast_numTriggerExtensions.addr = alloca <4 x i32> addrspace(1)*, align 4 ; <<4 x i32> addrspace(1)**> [#uses=2]
  %blast_numHits.addr = alloca <4 x i32> addrspace(1)*, align 4 ; <<4 x i32> addrspace(1)**> [#uses=2]
  %hitMatrix_furthestp.addr = alloca <4 x i32> addrspace(1)*, align 4 ; <<4 x i32> addrspace(1)**> [#uses=2]
  %hitMatrix_offsetp.addr = alloca <4 x i32> addrspace(1)*, align 4 ; <<4 x i32> addrspace(1)**> [#uses=2]
  %ungappedExtension_extensionsp.addr = alloca %struct.ungappedExtension addrspace(1)*, align 4 ; <%struct.ungappedExtension addrspace(1)**> [#uses=4]
  %nTotalSequenceNum.addr = alloca <4 x i32>, align 16 ; <<4 x i32>*> [#uses=2]
  %subject = alloca i8*, align 4                  ; <i8**> [#uses=4]
  %sequenceEnd = alloca i8*, align 4              ; <i8**> [#uses=3]
  %address = alloca i8*, align 4                  ; <i8**> [#uses=16]
  %subjectOffset = alloca <4 x i32>, align 16     ; <<4 x i32>*> [#uses=2]
  %count = alloca <4 x i32>, align 16             ; <<4 x i32>*> [#uses=4]
  %currentWord = alloca i8, align 1               ; <i8*> [#uses=4]
  %currentBlock = alloca i8*, align 4             ; <i8**> [#uses=4]
  %currentGroupFP = alloca %struct.coordinate*, align 4 ; <%struct.coordinate**> [#uses=10]
  %wordLookupDFA_AddiPositions = alloca <2 x i32>*, align 4 ; <<2 x i32>**> [#uses=2]
  %queryOffsets = alloca <2 x i32>*, align 4      ; <<2 x i32>**> [#uses=9]
  %queryOffset = alloca <2 x i32>, align 8        ; <<2 x i32>*> [#uses=3]
  %ungappedExtension_current = alloca %struct.ungappedExtension*, align 4 ; <%struct.ungappedExtension**> [#uses=2]
  %ungappedExtension_additional = alloca %struct.ungappedExtension*, align 4 ; <%struct.ungappedExtension**> [#uses=2]
  %diagonal = alloca <4 x i32>, align 16          ; <<4 x i32>*> [#uses=2]
  %lastHitFP = alloca <4 x i32>*, align 4         ; <<4 x i32>**> [#uses=5]
  %ungappedExtension_subjectEndReachedFP = alloca <4 x i32>, align 16 ; <<4 x i32>*> [#uses=2]
  %hitMatrix_Local = alloca <4 x i32>*, align 4   ; <<4 x i32>**> [#uses=2]
  %sequenceCount = alloca <4 x i32>, align 16     ; <<4 x i32>*> [#uses=8]
  %distance = alloca <4 x i32>, align 16          ; <<4 x i32>*> [#uses=3]
  %tid = alloca i32, align 4                      ; <i32*> [#uses=8]
  %numOfTriggerExtensions = alloca <4 x i32>, align 16 ; <<4 x i32>*> [#uses=3]
  %tmp225 = alloca <2 x i32>, align 8             ; <<2 x i32>*> [#uses=2]
  %tmp234 = alloca <4 x i32>, align 16            ; <<4 x i32>*> [#uses=9]
  %tmp235 = alloca <4 x i32>, align 16            ; <<4 x i32>*> [#uses=5]
  %tmp265 = alloca <4 x i32>, align 16            ; <<4 x i32>*> [#uses=2]
  %tmp333 = alloca %struct.PSSMatrixFP, align 16  ; <%struct.PSSMatrixFP*> [#uses=2]
  store %struct.PSSMatrixFP addrspace(1)* %PSSMatrixFP, %struct.PSSMatrixFP addrspace(1)** %PSSMatrixFP.addr
  store <2 x i32> addrspace(1)* %matrixBody, <2 x i32> addrspace(1)** %matrixBody.addr
  store %struct.sequenceDataFP addrspace(1)* %sequenceDataFP, %struct.sequenceDataFP addrspace(1)** %sequenceDataFP.addr
  store i8 addrspace(1)* %sequence, i8 addrspace(1)** %sequence.addr
  store %struct.parameters addrspace(1)* %parametersFP, %struct.parameters addrspace(1)** %parametersFP.addr
  store %struct.coordinate addrspace(1)* %wordLookupDFA_groupsFP, %struct.coordinate addrspace(1)** %wordLookupDFA_groupsFP.addr
  store i8 addrspace(1)* %wordLookupDFAFP, i8 addrspace(1)** %wordLookupDFAFP.addr
  store <4 x i32> addrspace(1)* %blast_numUngappedExtensions, <4 x i32> addrspace(1)** %blast_numUngappedExtensions.addr
  store <4 x i32> addrspace(1)* %blast_numTriggerExtensions, <4 x i32> addrspace(1)** %blast_numTriggerExtensions.addr
  store <4 x i32> addrspace(1)* %blast_numHits, <4 x i32> addrspace(1)** %blast_numHits.addr
  store <4 x i32> addrspace(1)* %hitMatrix_furthestp, <4 x i32> addrspace(1)** %hitMatrix_furthestp.addr
  store <4 x i32> addrspace(1)* %hitMatrix_offsetp, <4 x i32> addrspace(1)** %hitMatrix_offsetp.addr
  store %struct.ungappedExtension addrspace(1)* %ungappedExtension_extensionsp, %struct.ungappedExtension addrspace(1)** %ungappedExtension_extensionsp.addr
  store <4 x i32> %nTotalSequenceNum, <4 x i32>* %nTotalSequenceNum.addr
  %call = call i32 @get_local_id(i32 1) nounwind  ; <i32> [#uses=1]
  %call1 = call i32 @get_local_size(i32 1) nounwind ; <i32> [#uses=1]
  %tmp = mul i32 %call, %call1                    ; <i32> [#uses=1]
  %call2 = call i32 @get_group_id(i32 1) nounwind ; <i32> [#uses=1]
  %tmp3 = add i32 %tmp, %call2                    ; <i32> [#uses=1]
  store i32 %tmp3, i32* %tid
  store <4 x i32> zeroinitializer, <4 x i32>* %numOfTriggerExtensions
  %tmp4 = load <4 x i32> addrspace(1)** %hitMatrix_furthestp.addr ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp5 = load <4 x i32> addrspace(1)** %hitMatrix_offsetp.addr ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp6 = load i32* %tid                          ; <i32> [#uses=1]
  %arrayidx = getelementptr <4 x i32> addrspace(1)* %tmp5, i32 %tmp6 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp7 = load <4 x i32> addrspace(1)* %arrayidx  ; <<4 x i32>> [#uses=1]
  %tmp8 = extractelement <4 x i32> %tmp7, i32 0   ; <i32> [#uses=1]
  %tmp9 = getelementptr <4 x i32> addrspace(1)* %tmp4, i32 %tmp8 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp10 = load %struct.PSSMatrixFP addrspace(1)** %PSSMatrixFP.addr ; <%struct.PSSMatrixFP addrspace(1)*> [#uses=1]
  %structele = getelementptr inbounds %struct.PSSMatrixFP addrspace(1)* %tmp10, i32 0, i32 0 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp11 = load <4 x i32> addrspace(1)* %structele ; <<4 x i32>> [#uses=1]
  %tmp12 = extractelement <4 x i32> %tmp11, i32 0 ; <i32> [#uses=1]
  %tmp13 = getelementptr <4 x i32> addrspace(1)* %tmp9, i32 %tmp12 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %conv = ptrtoint <4 x i32> addrspace(1)* %tmp13 to i32 ; <i32> [#uses=1]
  %conv14 = inttoptr i32 %conv to <4 x i32>*      ; <<4 x i32>*> [#uses=1]
  store <4 x i32>* %conv14, <4 x i32>** %hitMatrix_Local
  %tmp15 = load %struct.ungappedExtension addrspace(1)** %ungappedExtension_extensionsp.addr ; <%struct.ungappedExtension addrspace(1)*> [#uses=1]
  %structele16 = getelementptr inbounds %struct.ungappedExtension addrspace(1)* %tmp15, i32 0, i32 0 ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %structele17 = getelementptr inbounds %struct.coordinate addrspace(1)* %structele16, i32 0, i32 1 ; <<4 x i32> addrspace(1)*> [#uses=1]
  store <4 x i32> zeroinitializer, <4 x i32> addrspace(1)* %structele17
  %tmp18 = load %struct.ungappedExtension addrspace(1)** %ungappedExtension_extensionsp.addr ; <%struct.ungappedExtension addrspace(1)*> [#uses=1]
  %tmp19 = load i32* %tid                         ; <i32> [#uses=1]
  %tmp20 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele21 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp20, i32 0, i32 21 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp22 = load <4 x i32> addrspace(1)* %structele21 ; <<4 x i32>> [#uses=1]
  %tmp23 = extractelement <4 x i32> %tmp22, i32 0 ; <i32> [#uses=1]
  %tmp24 = mul i32 %tmp19, %tmp23                 ; <i32> [#uses=1]
  %tmp25 = getelementptr %struct.ungappedExtension addrspace(1)* %tmp18, i32 %tmp24 ; <%struct.ungappedExtension addrspace(1)*> [#uses=1]
  %conv26 = ptrtoint %struct.ungappedExtension addrspace(1)* %tmp25 to i32 ; <i32> [#uses=1]
  %conv27 = inttoptr i32 %conv26 to %struct.ungappedExtension* ; <%struct.ungappedExtension*> [#uses=1]
  store %struct.ungappedExtension* %conv27, %struct.ungappedExtension** %ungappedExtension_current
  %tmp28 = load %struct.ungappedExtension addrspace(1)** %ungappedExtension_extensionsp.addr ; <%struct.ungappedExtension addrspace(1)*> [#uses=1]
  %tmp29 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele30 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp29, i32 0, i32 22 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp31 = load <4 x i32> addrspace(1)* %structele30 ; <<4 x i32>> [#uses=1]
  %tmp32 = extractelement <4 x i32> %tmp31, i32 0 ; <i32> [#uses=1]
  %tmp33 = getelementptr %struct.ungappedExtension addrspace(1)* %tmp28, i32 %tmp32 ; <%struct.ungappedExtension addrspace(1)*> [#uses=1]
  %conv34 = ptrtoint %struct.ungappedExtension addrspace(1)* %tmp33 to i32 ; <i32> [#uses=1]
  %conv35 = inttoptr i32 %conv34 to %struct.ungappedExtension* ; <%struct.ungappedExtension*> [#uses=1]
  store %struct.ungappedExtension* %conv35, %struct.ungappedExtension** %ungappedExtension_additional
  %tmp36 = load i8 addrspace(1)** %wordLookupDFAFP.addr ; <i8 addrspace(1)*> [#uses=1]
  %tmp37 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele38 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp37, i32 0, i32 17 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp39 = load <4 x i32> addrspace(1)* %structele38 ; <<4 x i32>> [#uses=1]
  %tmp40 = extractelement <4 x i32> %tmp39, i32 0 ; <i32> [#uses=1]
  %tmp41 = getelementptr i8 addrspace(1)* %tmp36, i32 %tmp40 ; <i8 addrspace(1)*> [#uses=1]
  %conv42 = ptrtoint i8 addrspace(1)* %tmp41 to i32 ; <i32> [#uses=1]
  %conv43 = inttoptr i32 %conv42 to <2 x i32>*    ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %conv43, <2 x i32>** %wordLookupDFA_AddiPositions
  %tmp44 = load %struct.PSSMatrixFP addrspace(1)** %PSSMatrixFP.addr ; <%struct.PSSMatrixFP addrspace(1)*> [#uses=1]
  %structele45 = getelementptr inbounds %struct.PSSMatrixFP addrspace(1)* %tmp44, i32 0, i32 4 ; <<2 x i32>* addrspace(1)*> [#uses=1]
  %tmp46 = load <2 x i32> addrspace(1)** %matrixBody.addr ; <<2 x i32> addrspace(1)*> [#uses=1]
  %conv47 = ptrtoint <2 x i32> addrspace(1)* %tmp46 to i32 ; <i32> [#uses=1]
  %tmp48 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele49 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp48, i32 0, i32 1 ; <i8 addrspace(1)*> [#uses=1]
  %tmp50 = load i8 addrspace(1)* %structele49     ; <i8> [#uses=1]
  %conv51 = zext i8 %tmp50 to i32                 ; <i32> [#uses=1]
  %tmp52 = add i32 %conv47, %conv51               ; <i32> [#uses=1]
  %conv53 = inttoptr i32 %tmp52 to <2 x i32>*     ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %conv53, <2 x i32>* addrspace(1)* %structele45
  %tmp54 = load i32* %tid                         ; <i32> [#uses=4]
  %conv55 = insertelement <4 x i32> undef, i32 %tmp54, i32 0 ; <<4 x i32>> [#uses=1]
  %conv56 = insertelement <4 x i32> %conv55, i32 %tmp54, i32 1 ; <<4 x i32>> [#uses=1]
  %conv57 = insertelement <4 x i32> %conv56, i32 %tmp54, i32 2 ; <<4 x i32>> [#uses=1]
  %conv58 = insertelement <4 x i32> %conv57, i32 %tmp54, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv58, <4 x i32>* %sequenceCount
  br label %while.cond

return:                                           ; preds = %while.exit
  ret void

while.cond:                                       ; preds = %if.end, %entry
  %tmp59 = load <4 x i32>* %sequenceCount         ; <<4 x i32>> [#uses=1]
  %tmp60 = extractelement <4 x i32> %tmp59, i32 0 ; <i32> [#uses=1]
  %tmp61 = load <4 x i32>* %nTotalSequenceNum.addr ; <<4 x i32>> [#uses=1]
  %tmp62 = extractelement <4 x i32> %tmp61, i32 0 ; <i32> [#uses=1]
  %cmp = icmp ult i32 %tmp60, %tmp62              ; <i1> [#uses=1]
  br i1 %cmp, label %while.body, label %while.exit

while.exit:                                       ; preds = %while.cond
  %tmp385 = load <4 x i32> addrspace(1)** %blast_numTriggerExtensions.addr ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp386 = load i32* %tid                        ; <i32> [#uses=1]
  %arrayidx387 = getelementptr <4 x i32> addrspace(1)* %tmp385, i32 %tmp386 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp388 = load <4 x i32>* %numOfTriggerExtensions ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp388, <4 x i32> addrspace(1)* %arrayidx387
  br label %return

while.body:                                       ; preds = %while.cond
  %tmp63 = load i8 addrspace(1)** %sequence.addr  ; <i8 addrspace(1)*> [#uses=1]
  %tmp64 = load %struct.sequenceDataFP addrspace(1)** %sequenceDataFP.addr ; <%struct.sequenceDataFP addrspace(1)*> [#uses=1]
  %tmp65 = load <4 x i32>* %sequenceCount         ; <<4 x i32>> [#uses=1]
  %tmp66 = extractelement <4 x i32> %tmp65, i32 0 ; <i32> [#uses=1]
  %arrayidx67 = getelementptr %struct.sequenceDataFP addrspace(1)* %tmp64, i32 %tmp66 ; <%struct.sequenceDataFP addrspace(1)*> [#uses=1]
  %structele68 = getelementptr inbounds %struct.sequenceDataFP addrspace(1)* %arrayidx67, i32 0, i32 4 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp69 = load <4 x i32> addrspace(1)* %structele68 ; <<4 x i32>> [#uses=1]
  %tmp70 = extractelement <4 x i32> %tmp69, i32 0 ; <i32> [#uses=1]
  %tmp71 = getelementptr i8 addrspace(1)* %tmp63, i32 %tmp70 ; <i8 addrspace(1)*> [#uses=1]
  %conv72 = ptrtoint i8 addrspace(1)* %tmp71 to i32 ; <i32> [#uses=1]
  %conv73 = inttoptr i32 %conv72 to i8*           ; <i8*> [#uses=2]
  store i8* %conv73, i8** %address
  store i8* %conv73, i8** %subject
  %tmp74 = load %struct.sequenceDataFP addrspace(1)** %sequenceDataFP.addr ; <%struct.sequenceDataFP addrspace(1)*> [#uses=1]
  %tmp75 = load <4 x i32>* %sequenceCount         ; <<4 x i32>> [#uses=1]
  %tmp76 = extractelement <4 x i32> %tmp75, i32 0 ; <i32> [#uses=1]
  %arrayidx77 = getelementptr %struct.sequenceDataFP addrspace(1)* %tmp74, i32 %tmp76 ; <%struct.sequenceDataFP addrspace(1)*> [#uses=1]
  %structele78 = getelementptr inbounds %struct.sequenceDataFP addrspace(1)* %arrayidx77, i32 0, i32 0 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp79 = load <4 x i32> addrspace(1)* %structele78 ; <<4 x i32>> [#uses=1]
  %tmp80 = extractelement <4 x i32> %tmp79, i32 0 ; <i32> [#uses=1]
  %tmp81 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele82 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp81, i32 0, i32 0 ; <i8 addrspace(1)*> [#uses=1]
  %tmp83 = load i8 addrspace(1)* %structele82     ; <i8> [#uses=1]
  %conv84 = sext i8 %tmp83 to i32                 ; <i32> [#uses=1]
  %cmp85 = icmp uge i32 %tmp80, %conv84           ; <i1> [#uses=1]
  br i1 %cmp85, label %if.then, label %if.end

if.end:                                           ; preds = %__T148605452, %while.body
  %tmp376 = load <4 x i32>* %sequenceCount        ; <<4 x i32>> [#uses=1]
  %call377 = call i32 @get_global_size(i32 1) nounwind ; <i32> [#uses=1]
  %call378 = call i32 @get_local_size(i32 1) nounwind ; <i32> [#uses=1]
  %tmp379 = mul i32 %call377, %call378            ; <i32> [#uses=4]
  %conv380 = insertelement <4 x i32> undef, i32 %tmp379, i32 0 ; <<4 x i32>> [#uses=1]
  %conv381 = insertelement <4 x i32> %conv380, i32 %tmp379, i32 1 ; <<4 x i32>> [#uses=1]
  %conv382 = insertelement <4 x i32> %conv381, i32 %tmp379, i32 2 ; <<4 x i32>> [#uses=1]
  %conv383 = insertelement <4 x i32> %conv382, i32 %tmp379, i32 3 ; <<4 x i32>> [#uses=1]
  %tmp384 = add <4 x i32> %tmp376, %conv383       ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp384, <4 x i32>* %sequenceCount
  br label %while.cond

if.then:                                          ; preds = %while.body
  %tmp86 = load %struct.coordinate addrspace(1)** %wordLookupDFA_groupsFP.addr ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %conv87 = ptrtoint %struct.coordinate addrspace(1)* %tmp86 to i32 ; <i32> [#uses=1]
  %conv88 = inttoptr i32 %conv87 to %struct.coordinate* ; <%struct.coordinate*> [#uses=1]
  store %struct.coordinate* %conv88, %struct.coordinate** %currentGroupFP
  store <4 x i32> <i32 1, i32 1, i32 1, i32 1>, <4 x i32>* %count
  br label %while.cond89

while.cond89:                                     ; preds = %if.end107, %if.then
  %tmp90 = load <4 x i32>* %count                 ; <<4 x i32>> [#uses=1]
  %tmp91 = extractelement <4 x i32> %tmp90, i32 0 ; <i32> [#uses=1]
  %tmp92 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele93 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp92, i32 0, i32 0 ; <i8 addrspace(1)*> [#uses=1]
  %tmp94 = load i8 addrspace(1)* %structele93     ; <i8> [#uses=1]
  %conv95 = sext i8 %tmp94 to i32                 ; <i32> [#uses=1]
  %cmp96 = icmp slt i32 %tmp91, %conv95           ; <i1> [#uses=1]
  br i1 %cmp96, label %while.body98, label %while.exit97

while.exit97:                                     ; preds = %while.cond89
  %tmp133 = load i8** %subject                    ; <i8*> [#uses=1]
  %tmp134 = load %struct.sequenceDataFP addrspace(1)** %sequenceDataFP.addr ; <%struct.sequenceDataFP addrspace(1)*> [#uses=1]
  %tmp135 = load <4 x i32>* %sequenceCount        ; <<4 x i32>> [#uses=1]
  %tmp136 = extractelement <4 x i32> %tmp135, i32 0 ; <i32> [#uses=1]
  %arrayidx137 = getelementptr %struct.sequenceDataFP addrspace(1)* %tmp134, i32 %tmp136 ; <%struct.sequenceDataFP addrspace(1)*> [#uses=1]
  %structele138 = getelementptr inbounds %struct.sequenceDataFP addrspace(1)* %arrayidx137, i32 0, i32 0 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp139 = load <4 x i32> addrspace(1)* %structele138 ; <<4 x i32>> [#uses=1]
  %tmp140 = extractelement <4 x i32> %tmp139, i32 0 ; <i32> [#uses=1]
  %tmp141 = getelementptr i8* %tmp133, i32 %tmp140 ; <i8*> [#uses=1]
  store i8* %tmp141, i8** %sequenceEnd
  br label %while.cond142

while.body98:                                     ; preds = %while.cond89
  %tmp99 = load i8** %address                     ; <i8*> [#uses=1]
  %tmp100 = load i8* %tmp99                       ; <i8> [#uses=1]
  %conv101 = zext i8 %tmp100 to i32               ; <i32> [#uses=1]
  %tmp102 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele103 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp102, i32 0, i32 16 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp104 = load <4 x i32> addrspace(1)* %structele103 ; <<4 x i32>> [#uses=1]
  %tmp105 = extractelement <4 x i32> %tmp104, i32 0 ; <i32> [#uses=1]
  %cmp106 = icmp slt i32 %conv101, %tmp105        ; <i1> [#uses=1]
  br i1 %cmp106, label %if.then108, label %if.else

if.end107:                                        ; preds = %if.else, %if.then108
  %tmp129 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp130 = getelementptr i8* %tmp129, i32 1      ; <i8*> [#uses=1]
  store i8* %tmp130, i8** %address
  %tmp131 = load <4 x i32>* %count                ; <<4 x i32>> [#uses=1]
  %tmp132 = add <4 x i32> %tmp131, <i32 1, i32 1, i32 1, i32 1> ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp132, <4 x i32>* %count
  br label %while.cond89

if.then108:                                       ; preds = %while.body98
  %tmp109 = load %struct.coordinate addrspace(1)** %wordLookupDFA_groupsFP.addr ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %tmp110 = load %struct.coordinate** %currentGroupFP ; <%struct.coordinate*> [#uses=1]
  %structele111 = getelementptr inbounds %struct.coordinate* %tmp110, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp112 = load <4 x i32>* %structele111         ; <<4 x i32>> [#uses=1]
  %tmp113 = extractelement <4 x i32> %tmp112, i32 0 ; <i32> [#uses=1]
  %tmp114 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp115 = load i8* %tmp114                      ; <i8> [#uses=1]
  %conv116 = zext i8 %tmp115 to i32               ; <i32> [#uses=1]
  %tmp117 = add i32 %tmp113, %conv116             ; <i32> [#uses=1]
  %tmp118 = getelementptr %struct.coordinate addrspace(1)* %tmp109, i32 %tmp117 ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %conv119 = ptrtoint %struct.coordinate addrspace(1)* %tmp118 to i32 ; <i32> [#uses=1]
  %conv120 = inttoptr i32 %conv119 to %struct.coordinate* ; <%struct.coordinate*> [#uses=1]
  store %struct.coordinate* %conv120, %struct.coordinate** %currentGroupFP
  br label %if.end107

if.else:                                          ; preds = %while.body98
  %tmp121 = load %struct.coordinate addrspace(1)** %wordLookupDFA_groupsFP.addr ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %tmp122 = load %struct.coordinate** %currentGroupFP ; <%struct.coordinate*> [#uses=1]
  %structele123 = getelementptr inbounds %struct.coordinate* %tmp122, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp124 = load <4 x i32>* %structele123         ; <<4 x i32>> [#uses=1]
  %tmp125 = extractelement <4 x i32> %tmp124, i32 0 ; <i32> [#uses=1]
  %tmp126 = getelementptr %struct.coordinate addrspace(1)* %tmp121, i32 %tmp125 ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %conv127 = ptrtoint %struct.coordinate addrspace(1)* %tmp126 to i32 ; <i32> [#uses=1]
  %conv128 = inttoptr i32 %conv127 to %struct.coordinate* ; <%struct.coordinate*> [#uses=1]
  store %struct.coordinate* %conv128, %struct.coordinate** %currentGroupFP
  br label %if.end107

while.cond142:                                    ; preds = %if.end201, %while.exit97
  %tmp143 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp144 = load i8** %sequenceEnd                ; <i8*> [#uses=1]
  %cmp145 = icmp ult i8* %tmp143, %tmp144         ; <i1> [#uses=1]
  br i1 %cmp145, label %while.body147, label %while.exit146

while.exit146:                                    ; preds = %while.cond142
  br label %__T148605452

while.body147:                                    ; preds = %while.cond142
  %tmp148 = load i8 addrspace(1)** %wordLookupDFAFP.addr ; <i8 addrspace(1)*> [#uses=1]
  %tmp149 = load %struct.coordinate** %currentGroupFP ; <%struct.coordinate*> [#uses=1]
  %structele150 = getelementptr inbounds %struct.coordinate* %tmp149, i32 0, i32 0 ; <<4 x i32>*> [#uses=1]
  %tmp151 = load <4 x i32>* %structele150         ; <<4 x i32>> [#uses=1]
  %tmp152 = extractelement <4 x i32> %tmp151, i32 0 ; <i32> [#uses=1]
  %tmp153 = getelementptr i8 addrspace(1)* %tmp148, i32 %tmp152 ; <i8 addrspace(1)*> [#uses=1]
  %conv154 = ptrtoint i8 addrspace(1)* %tmp153 to i32 ; <i32> [#uses=1]
  %conv155 = inttoptr i32 %conv154 to i8*         ; <i8*> [#uses=1]
  store i8* %conv155, i8** %currentBlock
  %tmp156 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp157 = load i8* %tmp156                      ; <i8> [#uses=1]
  %conv158 = zext i8 %tmp157 to i32               ; <i32> [#uses=1]
  %tmp159 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele160 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp159, i32 0, i32 16 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp161 = load <4 x i32> addrspace(1)* %structele160 ; <<4 x i32>> [#uses=1]
  %tmp162 = extractelement <4 x i32> %tmp161, i32 0 ; <i32> [#uses=1]
  %cmp163 = icmp slt i32 %conv158, %tmp162        ; <i1> [#uses=1]
  br i1 %cmp163, label %if.then165, label %if.else166

if.end164:                                        ; preds = %if.end188, %if.then165
  %tmp200 = load i8* %currentWord                 ; <i8> [#uses=1]
  %tobool = icmp ne i8 %tmp200, 0                 ; <i1> [#uses=1]
  br i1 %tobool, label %if.then202, label %if.end201

if.then165:                                       ; preds = %while.body147
  %tmp167 = load i8** %currentBlock               ; <i8*> [#uses=1]
  %tmp168 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp169 = load i8* %tmp168                      ; <i8> [#uses=1]
  %conv170 = zext i8 %tmp169 to i32               ; <i32> [#uses=1]
  %arrayidx171 = getelementptr i8* %tmp167, i32 %conv170 ; <i8*> [#uses=1]
  %tmp172 = load i8* %arrayidx171                 ; <i8> [#uses=1]
  store i8 %tmp172, i8* %currentWord
  %tmp173 = load %struct.coordinate addrspace(1)** %wordLookupDFA_groupsFP.addr ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %tmp174 = load %struct.coordinate** %currentGroupFP ; <%struct.coordinate*> [#uses=1]
  %structele175 = getelementptr inbounds %struct.coordinate* %tmp174, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp176 = load <4 x i32>* %structele175         ; <<4 x i32>> [#uses=1]
  %tmp177 = extractelement <4 x i32> %tmp176, i32 0 ; <i32> [#uses=1]
  %tmp178 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp179 = load i8* %tmp178                      ; <i8> [#uses=1]
  %conv180 = zext i8 %tmp179 to i32               ; <i32> [#uses=1]
  %tmp181 = add i32 %tmp177, %conv180             ; <i32> [#uses=1]
  %tmp182 = getelementptr %struct.coordinate addrspace(1)* %tmp173, i32 %tmp181 ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %conv183 = ptrtoint %struct.coordinate addrspace(1)* %tmp182 to i32 ; <i32> [#uses=1]
  %conv184 = inttoptr i32 %conv183 to %struct.coordinate* ; <%struct.coordinate*> [#uses=1]
  store %struct.coordinate* %conv184, %struct.coordinate** %currentGroupFP
  br label %if.end164

if.else166:                                       ; preds = %while.body147
  %tmp185 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp186 = load i8** %sequenceEnd                ; <i8*> [#uses=1]
  %cmp187 = icmp uge i8* %tmp185, %tmp186         ; <i1> [#uses=1]
  br i1 %cmp187, label %if.then189, label %if.end188

if.end188:                                        ; preds = %if.else166
  %tmp190 = load i8** %currentBlock               ; <i8*> [#uses=1]
  %tmp191 = load i8* %tmp190                      ; <i8> [#uses=1]
  store i8 %tmp191, i8* %currentWord
  %tmp192 = load %struct.coordinate addrspace(1)** %wordLookupDFA_groupsFP.addr ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %tmp193 = load %struct.coordinate** %currentGroupFP ; <%struct.coordinate*> [#uses=1]
  %structele194 = getelementptr inbounds %struct.coordinate* %tmp193, i32 0, i32 1 ; <<4 x i32>*> [#uses=1]
  %tmp195 = load <4 x i32>* %structele194         ; <<4 x i32>> [#uses=1]
  %tmp196 = extractelement <4 x i32> %tmp195, i32 0 ; <i32> [#uses=1]
  %tmp197 = getelementptr %struct.coordinate addrspace(1)* %tmp192, i32 %tmp196 ; <%struct.coordinate addrspace(1)*> [#uses=1]
  %conv198 = ptrtoint %struct.coordinate addrspace(1)* %tmp197 to i32 ; <i32> [#uses=1]
  %conv199 = inttoptr i32 %conv198 to %struct.coordinate* ; <%struct.coordinate*> [#uses=1]
  store %struct.coordinate* %conv199, %struct.coordinate** %currentGroupFP
  br label %if.end164

if.then189:                                       ; preds = %if.else166
  br label %__T148605452

__T148605452:                                     ; preds = %while.exit146, %if.then189
  br label %if.end

if.end201:                                        ; preds = %do.exit, %if.end164
  %tmp374 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp375 = getelementptr i8* %tmp374, i32 1      ; <i8*> [#uses=1]
  store i8* %tmp375, i8** %address
  br label %while.cond142

if.then202:                                       ; preds = %if.end164
  %tmp203 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp204 = load i8** %subject                    ; <i8*> [#uses=1]
  %sub.ptr.lhs.cast = ptrtoint i8* %tmp203 to i32 ; <i32> [#uses=1]
  %sub.ptr.rhs.cast = ptrtoint i8* %tmp204 to i32 ; <i32> [#uses=1]
  %sub.ptr.sub = sub i32 %sub.ptr.lhs.cast, %sub.ptr.rhs.cast ; <i32> [#uses=4]
  %conv205 = insertelement <4 x i32> undef, i32 %sub.ptr.sub, i32 0 ; <<4 x i32>> [#uses=1]
  %conv206 = insertelement <4 x i32> %conv205, i32 %sub.ptr.sub, i32 1 ; <<4 x i32>> [#uses=1]
  %conv207 = insertelement <4 x i32> %conv206, i32 %sub.ptr.sub, i32 2 ; <<4 x i32>> [#uses=1]
  %conv208 = insertelement <4 x i32> %conv207, i32 %sub.ptr.sub, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv208, <4 x i32>* %subjectOffset
  %tmp209 = load i8** %currentBlock               ; <i8*> [#uses=1]
  %conv210 = bitcast i8* %tmp209 to <2 x i32>*    ; <<2 x i32>*> [#uses=1]
  %tmp211 = load i8* %currentWord                 ; <i8> [#uses=1]
  %conv212 = zext i8 %tmp211 to i32               ; <i32> [#uses=1]
  %sub.ptr.neg = sub i32 0, %conv212              ; <i32> [#uses=1]
  %tmp213 = getelementptr <2 x i32>* %conv210, i32 %sub.ptr.neg ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp213, <2 x i32>** %queryOffsets
  %tmp214 = load <2 x i32>** %queryOffsets        ; <<2 x i32>*> [#uses=1]
  %tmp215 = load <2 x i32>* %tmp214               ; <<2 x i32>> [#uses=1]
  %tmp216 = extractelement <2 x i32> %tmp215, i32 0 ; <i32> [#uses=1]
  %tobool217 = icmp ne i32 %tmp216, 0             ; <i1> [#uses=1]
  %lnot = xor i1 %tobool217, true                 ; <i1> [#uses=1]
  br i1 %lnot, label %if.then219, label %if.end218

if.end218:                                        ; preds = %if.then219, %if.then202
  br label %do.body

if.then219:                                       ; preds = %if.then202
  %tmp220 = load <2 x i32>** %wordLookupDFA_AddiPositions ; <<2 x i32>*> [#uses=1]
  %tmp221 = load <2 x i32>** %queryOffsets        ; <<2 x i32>*> [#uses=1]
  %tmp222 = getelementptr <2 x i32>* %tmp221, i32 1 ; <<2 x i32>*> [#uses=1]
  %tmp223 = load <2 x i32>* %tmp222               ; <<2 x i32>> [#uses=1]
  %tmp224 = mul <2 x i32> %tmp223, <i32 32767, i32 32767> ; <<2 x i32>> [#uses=1]
  store <2 x i32> %tmp224, <2 x i32>* %tmp225
  %tmp226 = load <2 x i32>* %tmp225               ; <<2 x i32>> [#uses=1]
  %tmp227 = extractelement <2 x i32> %tmp226, i32 0 ; <i32> [#uses=1]
  %tmp228 = getelementptr <2 x i32>* %tmp220, i32 %tmp227 ; <<2 x i32>*> [#uses=1]
  %tmp229 = load <2 x i32>** %queryOffsets        ; <<2 x i32>*> [#uses=1]
  %tmp230 = getelementptr <2 x i32>* %tmp229, i32 2 ; <<2 x i32>*> [#uses=1]
  %tmp231 = load <2 x i32>* %tmp230               ; <<2 x i32>> [#uses=1]
  %tmp232 = extractelement <2 x i32> %tmp231, i32 0 ; <i32> [#uses=1]
  %tmp233 = getelementptr <2 x i32>* %tmp228, i32 %tmp232 ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp233, <2 x i32>** %queryOffsets
  br label %if.end218

do.exit:                                          ; preds = %do.cond
  br label %if.end201

do.body:                                          ; preds = %do.cond, %if.end218
  %tmp236 = load <2 x i32>** %queryOffsets        ; <<2 x i32>*> [#uses=1]
  %tmp237 = load <2 x i32>* %tmp236               ; <<2 x i32>> [#uses=1]
  store <2 x i32> %tmp237, <2 x i32>* %queryOffset
  %tmp238 = load <4 x i32>* %subjectOffset        ; <<4 x i32>> [#uses=1]
  %tmp239 = load <2 x i32>* %queryOffset          ; <<2 x i32>> [#uses=1]
  %tmp240 = extractelement <2 x i32> %tmp239, i32 0 ; <i32> [#uses=4]
  %conv241 = insertelement <4 x i32> undef, i32 %tmp240, i32 0 ; <<4 x i32>> [#uses=1]
  %conv242 = insertelement <4 x i32> %conv241, i32 %tmp240, i32 1 ; <<4 x i32>> [#uses=1]
  %conv243 = insertelement <4 x i32> %conv242, i32 %tmp240, i32 2 ; <<4 x i32>> [#uses=1]
  %conv244 = insertelement <4 x i32> %conv243, i32 %tmp240, i32 3 ; <<4 x i32>> [#uses=1]
  %tmp245 = sub <4 x i32> %tmp238, %conv244       ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp245, <4 x i32>* %diagonal
  %tmp246 = load <4 x i32>** %hitMatrix_Local     ; <<4 x i32>*> [#uses=1]
  %tmp247 = load <4 x i32>* %diagonal             ; <<4 x i32>> [#uses=1]
  %tmp248 = extractelement <4 x i32> %tmp247, i32 0 ; <i32> [#uses=1]
  %tmp249 = getelementptr <4 x i32>* %tmp246, i32 %tmp248 ; <<4 x i32>*> [#uses=1]
  store <4 x i32>* %tmp249, <4 x i32>** %lastHitFP
  %tmp250 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp251 = load i8 addrspace(1)** %sequence.addr ; <i8 addrspace(1)*> [#uses=1]
  %conv252 = ptrtoint i8 addrspace(1)* %tmp251 to i32 ; <i32> [#uses=1]
  %conv253 = inttoptr i32 %conv252 to i8*         ; <i8*> [#uses=1]
  %sub.ptr.lhs.cast254 = ptrtoint i8* %tmp250 to i32 ; <i32> [#uses=1]
  %sub.ptr.rhs.cast255 = ptrtoint i8* %conv253 to i32 ; <i32> [#uses=1]
  %sub.ptr.sub256 = sub i32 %sub.ptr.lhs.cast254, %sub.ptr.rhs.cast255 ; <i32> [#uses=4]
  %conv257 = insertelement <4 x i32> undef, i32 %sub.ptr.sub256, i32 0 ; <<4 x i32>> [#uses=1]
  %conv258 = insertelement <4 x i32> %conv257, i32 %sub.ptr.sub256, i32 1 ; <<4 x i32>> [#uses=1]
  %conv259 = insertelement <4 x i32> %conv258, i32 %sub.ptr.sub256, i32 2 ; <<4 x i32>> [#uses=1]
  %conv260 = insertelement <4 x i32> %conv259, i32 %sub.ptr.sub256, i32 3 ; <<4 x i32>> [#uses=1]
  %tmp261 = load <4 x i32>** %lastHitFP           ; <<4 x i32>*> [#uses=1]
  %tmp262 = load <4 x i32>* %tmp261               ; <<4 x i32>> [#uses=1]
  %tmp263 = sub <4 x i32> %conv260, %tmp262       ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp263, <4 x i32>* %tmp235
  %tmp264 = load <4 x i32>* %tmp235               ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp264, <4 x i32>* %tmp265
  %tmp266 = load <4 x i32>* %tmp265               ; <<4 x i32>> [#uses=1]
  %tmp267 = extractelement <4 x i32> %tmp266, i32 0 ; <i32> [#uses=1]
  %tmp268 = load <4 x i32>* %tmp234               ; <<4 x i32>> [#uses=1]
  %tmp269 = insertelement <4 x i32> %tmp268, i32 %tmp267, i32 0 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp269, <4 x i32>* %tmp234
  %tmp270 = load <4 x i32>* %tmp235               ; <<4 x i32>> [#uses=1]
  %tmp271 = extractelement <4 x i32> %tmp270, i32 1 ; <i32> [#uses=1]
  %tmp272 = load <4 x i32>* %tmp234               ; <<4 x i32>> [#uses=1]
  %tmp273 = insertelement <4 x i32> %tmp272, i32 %tmp271, i32 1 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp273, <4 x i32>* %tmp234
  %tmp274 = load <4 x i32>* %tmp235               ; <<4 x i32>> [#uses=1]
  %tmp275 = extractelement <4 x i32> %tmp274, i32 2 ; <i32> [#uses=1]
  %tmp276 = load <4 x i32>* %tmp234               ; <<4 x i32>> [#uses=1]
  %tmp277 = insertelement <4 x i32> %tmp276, i32 %tmp275, i32 2 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp277, <4 x i32>* %tmp234
  %tmp278 = load <4 x i32>* %tmp235               ; <<4 x i32>> [#uses=1]
  %tmp279 = extractelement <4 x i32> %tmp278, i32 3 ; <i32> [#uses=1]
  %tmp280 = load <4 x i32>* %tmp234               ; <<4 x i32>> [#uses=1]
  %tmp281 = insertelement <4 x i32> %tmp280, i32 %tmp279, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp281, <4 x i32>* %tmp234
  %tmp282 = load <4 x i32>* %tmp234               ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp282, <4 x i32>* %distance
  %tmp283 = load <4 x i32>* %distance             ; <<4 x i32>> [#uses=1]
  %tmp284 = extractelement <4 x i32> %tmp283, i32 0 ; <i32> [#uses=1]
  %tmp285 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele286 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp285, i32 0, i32 20 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp287 = load <4 x i32> addrspace(1)* %structele286 ; <<4 x i32>> [#uses=1]
  %tmp288 = extractelement <4 x i32> %tmp287, i32 0 ; <i32> [#uses=1]
  %cmp289 = icmp sge i32 %tmp284, %tmp288         ; <i1> [#uses=1]
  br i1 %cmp289, label %if.then291, label %if.else292

do.cond:                                          ; preds = %if.end290
  %tmp370 = load <2 x i32>** %queryOffsets        ; <<2 x i32>*> [#uses=1]
  %tmp371 = load <2 x i32>* %tmp370               ; <<2 x i32>> [#uses=1]
  %tmp372 = extractelement <2 x i32> %tmp371, i32 0 ; <i32> [#uses=1]
  %tobool373 = icmp ne i32 %tmp372, 0             ; <i1> [#uses=1]
  br i1 %tobool373, label %do.body, label %do.exit

if.end290:                                        ; preds = %if.end312, %if.then291
  %tmp363 = load <2 x i32>** %queryOffsets        ; <<2 x i32>*> [#uses=1]
  %tmp364 = getelementptr <2 x i32>* %tmp363, i32 1 ; <<2 x i32>*> [#uses=1]
  store <2 x i32>* %tmp364, <2 x i32>** %queryOffsets
  %tmp365 = load <4 x i32> addrspace(1)** %blast_numHits.addr ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp366 = load i32* %tid                        ; <i32> [#uses=1]
  %arrayidx367 = getelementptr <4 x i32> addrspace(1)* %tmp365, i32 %tmp366 ; <<4 x i32> addrspace(1)*> [#uses=2]
  %tmp368 = load <4 x i32> addrspace(1)* %arrayidx367 ; <<4 x i32>> [#uses=1]
  %tmp369 = add <4 x i32> %tmp368, <i32 1, i32 1, i32 1, i32 1> ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp369, <4 x i32> addrspace(1)* %arrayidx367
  br label %do.cond

if.then291:                                       ; preds = %do.body
  %tmp293 = load <4 x i32>** %lastHitFP           ; <<4 x i32>*> [#uses=1]
  %tmp294 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp295 = load i8 addrspace(1)** %sequence.addr ; <i8 addrspace(1)*> [#uses=1]
  %conv296 = ptrtoint i8 addrspace(1)* %tmp295 to i32 ; <i32> [#uses=1]
  %conv297 = inttoptr i32 %conv296 to i8*         ; <i8*> [#uses=1]
  %sub.ptr.lhs.cast298 = ptrtoint i8* %tmp294 to i32 ; <i32> [#uses=1]
  %sub.ptr.rhs.cast299 = ptrtoint i8* %conv297 to i32 ; <i32> [#uses=1]
  %sub.ptr.sub300 = sub i32 %sub.ptr.lhs.cast298, %sub.ptr.rhs.cast299 ; <i32> [#uses=4]
  %conv301 = insertelement <4 x i32> undef, i32 %sub.ptr.sub300, i32 0 ; <<4 x i32>> [#uses=1]
  %conv302 = insertelement <4 x i32> %conv301, i32 %sub.ptr.sub300, i32 1 ; <<4 x i32>> [#uses=1]
  %conv303 = insertelement <4 x i32> %conv302, i32 %sub.ptr.sub300, i32 2 ; <<4 x i32>> [#uses=1]
  %conv304 = insertelement <4 x i32> %conv303, i32 %sub.ptr.sub300, i32 3 ; <<4 x i32>> [#uses=1]
  store <4 x i32> %conv304, <4 x i32>* %tmp293
  br label %if.end290

if.else292:                                       ; preds = %do.body
  %tmp305 = load <4 x i32>* %distance             ; <<4 x i32>> [#uses=1]
  %tmp306 = extractelement <4 x i32> %tmp305, i32 0 ; <i32> [#uses=1]
  %tmp307 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele308 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp307, i32 0, i32 2 ; <i8 addrspace(1)*> [#uses=1]
  %tmp309 = load i8 addrspace(1)* %structele308   ; <i8> [#uses=1]
  %conv310 = sext i8 %tmp309 to i32               ; <i32> [#uses=1]
  %cmp311 = icmp sge i32 %tmp306, %conv310        ; <i1> [#uses=1]
  br i1 %cmp311, label %if.then313, label %if.end312

if.end312:                                        ; preds = %if.then313, %if.else292
  br label %if.end290

if.then313:                                       ; preds = %if.else292
  %tmp314 = load <4 x i32> addrspace(1)** %blast_numUngappedExtensions.addr ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp315 = load i32* %tid                        ; <i32> [#uses=1]
  %arrayidx316 = getelementptr <4 x i32> addrspace(1)* %tmp314, i32 %tmp315 ; <<4 x i32> addrspace(1)*> [#uses=2]
  %tmp317 = load <4 x i32> addrspace(1)* %arrayidx316 ; <<4 x i32>> [#uses=1]
  %tmp318 = add <4 x i32> %tmp317, <i32 1, i32 1, i32 1, i32 1> ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp318, <4 x i32> addrspace(1)* %arrayidx316
  %tmp319 = load i8 addrspace(1)** %sequence.addr ; <i8 addrspace(1)*> [#uses=1]
  %conv320 = ptrtoint i8 addrspace(1)* %tmp319 to i32 ; <i32> [#uses=1]
  %conv321 = inttoptr i32 %conv320 to i8*         ; <i8*> [#uses=1]
  %tmp322 = load <2 x i32>* %queryOffset          ; <<2 x i32>> [#uses=1]
  %tmp323 = extractelement <2 x i32> %tmp322, i32 0 ; <i32> [#uses=4]
  %conv324 = insertelement <4 x i32> undef, i32 %tmp323, i32 0 ; <<4 x i32>> [#uses=1]
  %conv325 = insertelement <4 x i32> %conv324, i32 %tmp323, i32 1 ; <<4 x i32>> [#uses=1]
  %conv326 = insertelement <4 x i32> %conv325, i32 %tmp323, i32 2 ; <<4 x i32>> [#uses=1]
  %conv327 = insertelement <4 x i32> %conv326, i32 %tmp323, i32 3 ; <<4 x i32>> [#uses=1]
  %tmp328 = load i8** %address                    ; <i8*> [#uses=1]
  %tmp329 = load <4 x i32>** %lastHitFP           ; <<4 x i32>*> [#uses=1]
  %tmp330 = load <4 x i32>* %tmp329               ; <<4 x i32>> [#uses=1]
  %tmp331 = load %struct.PSSMatrixFP addrspace(1)** %PSSMatrixFP.addr ; <%struct.PSSMatrixFP addrspace(1)*> [#uses=1]
  %tmp332 = load %struct.PSSMatrixFP addrspace(1)* %tmp331 ; <%struct.PSSMatrixFP> [#uses=1]
  store %struct.PSSMatrixFP %tmp332, %struct.PSSMatrixFP* %tmp333
  %tmp334 = load i8** %subject                    ; <i8*> [#uses=1]
  %tmp335 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele336 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp335, i32 0, i32 1 ; <i8 addrspace(1)*> [#uses=1]
  %tmp337 = load i8 addrspace(1)* %structele336   ; <i8> [#uses=1]
  %tmp338 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele339 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp338, i32 0, i32 18 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp340 = load <4 x i32> addrspace(1)* %structele339 ; <<4 x i32>> [#uses=1]
  %tmp341 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele342 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp341, i32 0, i32 19 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp343 = load <4 x i32> addrspace(1)* %structele342 ; <<4 x i32>> [#uses=1]
  %tmp344 = load %struct.parameters addrspace(1)** %parametersFP.addr ; <%struct.parameters addrspace(1)*> [#uses=1]
  %structele345 = getelementptr inbounds %struct.parameters addrspace(1)* %tmp344, i32 0, i32 21 ; <<4 x i32> addrspace(1)*> [#uses=1]
  %tmp346 = load <4 x i32> addrspace(1)* %structele345 ; <<4 x i32>> [#uses=1]
  %tmp347 = extractelement <4 x i32> %tmp346, i32 0 ; <i32> [#uses=4]
  %conv348 = insertelement <4 x i32> undef, i32 %tmp347, i32 0 ; <<4 x i32>> [#uses=1]
  %conv349 = insertelement <4 x i32> %conv348, i32 %tmp347, i32 1 ; <<4 x i32>> [#uses=1]
  %conv350 = insertelement <4 x i32> %conv349, i32 %tmp347, i32 2 ; <<4 x i32>> [#uses=1]
  %conv351 = insertelement <4 x i32> %conv350, i32 %tmp347, i32 3 ; <<4 x i32>> [#uses=1]
  %tmp352 = load %struct.ungappedExtension** %ungappedExtension_current ; <%struct.ungappedExtension*> [#uses=1]
  %tmp353 = load %struct.ungappedExtension** %ungappedExtension_additional ; <%struct.ungappedExtension*> [#uses=1]
  %tmp354 = load <4 x i32>* %sequenceCount        ; <<4 x i32>> [#uses=1]
  %tmp355 = load i32* %tid                        ; <i32> [#uses=4]
  %conv356 = insertelement <4 x i32> undef, i32 %tmp355, i32 0 ; <<4 x i32>> [#uses=1]
  %conv357 = insertelement <4 x i32> %conv356, i32 %tmp355, i32 1 ; <<4 x i32>> [#uses=1]
  %conv358 = insertelement <4 x i32> %conv357, i32 %tmp355, i32 2 ; <<4 x i32>> [#uses=1]
  %conv359 = insertelement <4 x i32> %conv358, i32 %tmp355, i32 3 ; <<4 x i32>> [#uses=1]
  %call360 = call %struct.ungappedExtension addrspace(1)* @ungappedExtension_twoHitExtendD(i8* %conv321, <4 x i32> %conv327, i8* %tmp328, <4 x i32> %tmp330, %struct.PSSMatrixFP* byval %tmp333, i8* %tmp334, <4 x i32>* %ungappedExtension_subjectEndReachedFP, i8 zeroext %tmp337, <4 x i32> %tmp340, <4 x i32> %tmp343, <4 x i32> %conv351, %struct.ungappedExtension* %tmp352, %struct.ungappedExtension* %tmp353, <4 x i32>* %numOfTriggerExtensions, <4 x i32> %tmp354, <4 x i32> %conv359) nounwind ; <%struct.ungappedExtension addrspace(1)*> [#uses=0]
  %tmp361 = load <4 x i32>** %lastHitFP           ; <<4 x i32>*> [#uses=1]
  %tmp362 = load <4 x i32>* %ungappedExtension_subjectEndReachedFP ; <<4 x i32>> [#uses=1]
  store <4 x i32> %tmp362, <4 x i32>* %tmp361
  br label %if.end312
}
