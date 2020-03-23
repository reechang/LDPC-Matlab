    function codeWord = ldpcEncoding(H,u)
    %% This program allow the encoding of an LDPC code using the parity-check matrix even if this matrix is not 
    %% a full rank matrix. 
    
    %% Ref: The paper entitled "Modified Approximate Lower Triangular Encoding of LDPC Codes," 2015 International Conference
    %% on Advances in Computer Engineering and Applications (ICACEA) IMS Engineering College, Ghaziabad, India (IEEE).
    
    %% Input variables:
    %% H:      Parity-check matrix.
    %% u:      Information bit vector. 
    
    %% Output variable:
    %% codeWord:      The produced codeword for the information bit vector u.
    
    
    
    %%*******Remove all zero rows********  
    H = H(any(H,2),:);
    %%***********************************
    
    Horg = H;%%Original matrix H   
    %%**************The initial step for the algorithm  H = [A T] ***********************
    
    %%******Check the sub matrix T ********
    m = size(H,1);
    n = size(H,2);
    
    A = [];
    T = [];
    swap = [];
    
    idx1 = 0;
    for j = n-m+1:n
      if (H(1,j)==1 & j==n-m+1)
        idx1 = 1;
        break; 
      elseif (H(1,j)==1 & j>n-m+1)
            idx1 = 1;
            H(:,[n-m+1,j]) = H(:,[j,n-m+1]); 
            swap = [swap; n-m+1 j];          
            break;
      endif
    endfor  
    %%*****************************************
    %%******Check the sub matrix A ************
    for j = 1:n-m
          if (H(1,j)==1 & idx1 == 0) 
            H(:,[j,n-m+1]) = H(:,[n-m+1,j]); 
            idx1 = 1;
            swap = [swap; j n-m+1];          
            break; 
        endif
    endfor 
    %%******************************************
    
    %%*************************************************************************
    
    %%***************************Beginning of the algorithm******************************
    for i = 1:m
      j = n-m+i;
        if i==1
          idxOnes = find(H(:,j)==1);
          for iOnes = 1:length(idxOnes)%%Addition of rows of the matrix that have 1 in a specific column  
            if (idxOnes(iOnes)>=i+1)
              H(idxOnes(iOnes),:) = mod(H(idxOnes(iOnes),:) + H(i,:),2);
            endif
          endfor
          %%*****Move the rows with all zeros to the end of H*****%%
             idx_zeros = find(all(H==0,2)); %Index of all zero rows.
             if(length(idx_zeros)>=1)
               H(idx_zeros,:)=[];%Remove those rows.%%This instruction don't work if H is in galois field
               H = [H; zeros(length(idx_zeros),n)];%Add all zero rows to the end of H.
             endif
          %%******************************************************%%       
         else         
      %%******Check the sub matrix T ********
        idx2 = 0;
        for j2 = n-m+i:n
          if (H(i,j2)==1 & j2==n-m+i)
            idx2 = 1;
            break; 
          elseif (H(i,j2)==1 & j2>n-m+i)
                idx2 = 1;
                H(:,[n-m+i,j2]) = H(:,[j2,n-m+i]); 
                swap = [swap; n-m+i j2];          
                break;
          endif
        endfor  
      %%*****************************************
      %%******Check the sub matrix A ************
        for j2 = 1:n-m
              if (H(i,j2)==1 & idx2 == 0) 
                H(:,[j2,n-m+i]) = H(:,[n-m+i,j2]); 
                idx2 = 1;
                swap = [swap; j2 n-m+i];         
                break; 
              endif
        endfor 
     %*****************************************
        idxOnes = find(H(:,j)==1);
        
          for iOnes = 1:length(idxOnes)%%Addition of rows of the matrix that have 1 in a specific column 
            if (idxOnes(iOnes)>=i+1)
              H(idxOnes(iOnes),:) = mod(H(idxOnes(iOnes),:) + H(i,:),2);
            endif
          endfor
         %%*********Move the all zero rows to the end of H********%%
             idx_zeros = find(all(H==0,2)); %Index of all zero rows.
             if(length(idx_zeros)>=1)
               H(idx_zeros,:)=[];%Remove those rows.
               H = [H; zeros(length(idx_zeros),n)];%Add all zero rows to the end of H.
             endif
          %%******************************************************%%       
         
      endif    
  endfor
  %%***********************************************************************************
  
  %%***Remove rows with all zeros from H and specify the matricies A and T*************
  idx_zeros = find(all(H==0,2)); %Index of all zero rows.
  r_swap = length(idx_zeros); 
    if r_swap>=1
      H(idx_zeros,:)=[];%Remove those rows.
      A  = [H(:,1:n-m) H(:,n-r_swap+1:n)]; 
      T = H(:,n-m+1:n-r_swap);         
    else
      H(idx_zeros,:)=[];%Remove those rows.
      A = H(:,1:n-m);
      T = H(:,n-m+1:n);      
    endif
        
  %%***********************************************************************************
  
  %%*******************Encoding process using backward substitution********************
    a = T;
    b = mod(A * u',2); 
    %%*********************BackwardSub algorithm: is not my algorithm**************************
    l = length(b);
    y(l,1) = mod(b(l)/a(l,l),2);
    for i = l-1:-1:1
      y(i,1)=mod((b(i)-a(i,i+1:l)*y(i+1:l,1))./a(i,i),2);%  Parity bit vector.
    end
    %%***********************************************************************************
    
    c = [u y']; % The codeword
    
    %%*****Apply the inverse permutation to the codeword(using swap and r)******
    if(r_swap>=1)
      c1 = c([1:n-m n-m+r_swap+1:n]);
      c2 = c(n-m+1:n-m+r_swap);
      c  = [c1 c2];
    endif
    if (~isempty(swap))
      for i=size(swap,1):-1:1
        c(:,[swap(i,1),swap(i,2)]) = c(:,[swap(i,2),swap(i,1)]);
      endfor      
    endif   
    %%**************************************************************************
    
  %%***********************************************************************************
    codeWord = c;
    endfunction
            
