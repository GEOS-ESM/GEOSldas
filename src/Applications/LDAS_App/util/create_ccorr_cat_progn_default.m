
% create string array to initialize ccorr_cat_progn in default namelist file
%
% reichle, 1 Jun 2005
% reichle, 6 Dec 2013 - updated for new "progn_pert_type"

sm={
'catdef' 
'rzexc'  
'srfexc' 
'snow'
'tc'
'ght(1)' 
'ght(2)' 
'ght(3)' 
'ght(4)' 
'ght(5)' 
'ght(6)' }

k=0;

for i=1:length(sm)
  for j=(i+1):length(sm)
    
    k=k+1;
    
    tmpstr = [ 'ccorr_progn_pert%', sm{i}, '%', sm{j}];
    
    os{k} = tmpstr;
    
    es{k} = ' = 0.';
    
    if (j==25)
      
      k=k+1;
      
      os{k} = '';
      es{k} = '';
      
    end
      
  end
end

sa = [ char(os) char(es) ];

diary tmp.txt

disp(sa)

diary off

% ========= EOF ====================================




