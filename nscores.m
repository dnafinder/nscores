function STATS=nscores(x,varargin)
% NSCORES Executes the Van der Waerden version of non parametric tests (Normal
% scores tests).
% Named for the Dutch mathematician Bartel Leendert van der Waerden, the Van der
% Waerden test is a statistical test that k population distribution functions
% are equal. The Van Der Waerden test converts the ranks to quantiles of the
% standard normal distribution. These are called normal scores and the test is
% computed from these normal scores. The standard ANOVA assumes that the errors
% (i.e., residuals) are normally distributed. If this normality assumption is
% not valid, an alternative is to use a non-parametric test. The advantage of
% the Van Der Waerden test is that it provides the high efficiency of the
% standard ANOVA analysis when the normality assumptions are in fact satisfied,
% but it also provides the robustness of the non-parametric test when the
% normality assumptions are not satisfied.
% This function compute the Normal Scores of 5 tests:
% - Levene, Mann-Whitney-Wilcoxon and Wilcoxon tests when there are 2 groups;
% - Kruskal-Wallis and Friedman test whene there are more than 2 groups.
%
% The function will use a GUI to select the proper test.
% Moreover, the GUI will ask which version of Normal score do you want to
% use: Blom, Tukey, Rankit, Van der Waerden
%
% Syntax: 	STATS=nscores(x,alpha)
%      
%     Inputs:
%           X - data matrix (Size of matrix must be n-by-2; data=column 1, group=column 2).
%           ALPHA - significance level (default = 0.05).
%
%     Outputs:
%           - Z
%           - p-value
%
%        If STATS nargout was specified the results will be stored in the STATS
%        struct.
%
%      Example: 
%
%                             Sample
%                   -------------------------
%                      1        2        3    
%                   -------------------------
%                      8.84     8.65     7.89
%                      9.92    10.7 	 9.16
%                      7.2     10.24	 7.34
%                      9.25     8.62	10.28
%                      9.45     9.94	 9.12
%                      9.14    10.55	 9.24
%                      9.99    10.13	 8.4
%                      9.21	    9.78	 8.6
%                      9.06	    9.01	 8.04
%                                		 8.45
%                                        9.51
%                                        8.15
%                                        7.69
%                   -------------------------
%
%       Data matrix must be:
% x=[8.84 9.92 7.20 9.25 9.45 9.14 9.99 9.21 9.06 8.65 10.70 10.24 8.62 9.94...
% 10.55 10.13 9.78 9.01 7.79 9.16 7.64 10.28 9.12 9.24 8.40 8.60 8.04 8.45...
% 9.51 8.15 7.69; repmat(1,1,9) repmat(2,1,9) repmat(3,1,13)]';
%
%           Calling on Matlab:
%           - nscores(x(1:18,:)) a GUI will ask you what test to perform among Levene, Mann-Whitney-Wilcoxon or Wilcoxon tests
%           - nscores(x(1:17,:)) a GUI will ask you which test to perform between Levene and Mann-Whitney-Wilcoxon tests
%           - nscores(x(1:27,:)) a GUI will ask you which test to perform between Kruskal-Wallis and Friedman tests
%           - nscores(x) a GUI will inform you that the only test that is possible to perform is the Kruskal-Wallis test
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2010). NSCORES: Normal scores version of several non-parametric tests.
% http://www.mathworks.com/matlabcentral/fileexchange/26855

%Input Error handling
p = inputParser;
addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','nonempty','ncols',2}));
addOptional(p,'alpha',0.05, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));
parse(p,x,varargin{:});
assert(all(x(:,2) == fix(x(:,2))),'Warning: all elements of column 2 of input matrix must be whole numbers')
assert(max(x(:,2))>1,'Warning: almost two groups are required...')
alpha=p.Results.alpha;
clear p

method={'Blom','Tukey','Rankit','Van der Waerden'};
type=listdlg('PromptString','Select a version:','ListSize',[300 150],...
    'Name','Disposable methods', 'SelectionMode','single',...
                      'ListString',method);

switch method{type}
    case 'Blom'
        c=3/8;
    case 'Tukey'
        c=1/3;
    case 'Rankit'
        c=1/2;
    case 'Van der Waerden'
        c=0;
end

if nargout
    STATS.method=method{type};
end

k=max(x(:,2));
n=zeros(1,k);
for I=1:k
    n(I)=length(x(x(:,2)==I));
end

if k==2
    if isequal(n./n(1),ones(size(n)))
        avtest={'Levene','Mann-Whitney-Wilcoxon (unpaired samples)','Wilcoxon test (paired samples)'};
        tst=listdlg('PromptString','Select the test that you want to perform:','ListSize',[300 150],...
            'Name','Disposable test', 'SelectionMode','single',...
            'ListString',avtest);
    else
        avtest={'Levene','Mann-Whitney-Wilcoxon (unpaired samples)'};
        tst=listdlg('PromptString','Select the test that you want to perform:','ListSize',[300 150],...
            'Name','Disposable test', 'SelectionMode','single',...
            'ListString',avtest);
    end
else
    if isequal(n./n(1),ones(size(n)))
        avtest={'Kruskal-Wallis (one way ANOVA)','Friedman test (two way ANOVA; without replication)'};
        tst=listdlg('PromptString','Select the test that you want to perform:','ListSize',[300 150],...
            'Name','Disposable test', 'SelectionMode','single',...
            'ListString',avtest);
    else
        avtest={'Kruskal-Wallis (one way ANOVA)'};
        tst=listdlg('PromptString','Select the test that you want to perform:','ListSize',[300 150],...
            'Name','Disposable test', 'SelectionMode','single',...
            'ListString',avtest);
    end
end

if nargout
    STATS.test=avtest{tst};
end

switch avtest{tst}
    case 'Levene' %Levene Test
        L=length(x); Matrix=zeros(L,8); %set the basic parameter
        Matrix(:,[1 2])=x;
        for I=1:k
            Matrix(Matrix(:,2)==I,3)=Matrix(Matrix(:,2)==I,1)-mean(Matrix(Matrix(:,2)==I,1));
        end
        Matrix(:,4)=tiedrank(Matrix(:,3)); %ranks
        Matrix(:,5)=(Matrix(:,4)-c)./(L+1-2*c); %transform ranks in quantile
        Matrix(:,6)=norminv(Matrix(:,5)); %normal scores
        Matrix(:,7)=Matrix(:,6).^2;
        Matrix(:,8)=Matrix(:,6).^4;
        [m,I]=min(n);
        A=[sum(Matrix(Matrix(:,2)==I,7)) sum(Matrix(:,7:8))];
        Z=(A(1)-(m/L*A(2)))/realsqrt(prod(n)/L/(L-1)*(A(3)-1/L*A(2)^2)); %Statistics
        clear Matrix L I m A 
    case 'Mann-Whitney-Wilcoxon (unpaired samples)' 
        L=length(x); Matrix=zeros(L,6); %set the basic parameter
        Matrix(:,[1 2])=x;
        Matrix(:,3)=tiedrank(Matrix(:,1)); %ranks
        Matrix(:,4)=(Matrix(:,3)-c)./(L+1-2*c); %transform ranks in quantile
        Matrix(:,5)=norminv(Matrix(:,4)); %normal scores
        Matrix(:,6)=Matrix(:,5).^2; %square of normal scores
        Zbar=zeros(1,k); %preallocate mean normal scores array
        for I=1:k
            Zbar(I)=mean(Matrix(Matrix(:,2)==I,5)); %mean normal scores
        end
        s2=sum(Matrix(:,6))/(L-1);
        Z=realsqrt(sum(n.*Zbar.^2)/s2); %Statistics
        clear Matrix L I Zbar s2
     case 'Kruskal-Wallis (one way ANOVA)' 
        L=length(x); Matrix=zeros(L,6); %set the basic parameter
        Matrix(:,[1 2])=x;
        Matrix(:,3)=tiedrank(Matrix(:,1)); %ranks
        Matrix(:,4)=(Matrix(:,3)-c)./(L+1-2*c); %transform ranks in quantile
        Matrix(:,5)=norminv(Matrix(:,4)); %normal scores
        Matrix(:,6)=Matrix(:,5).^2; %square of normal scores
        Zbar=zeros(1,k); %preallocate mean normal scores array
        for I=1:k
            Zbar(I)=mean(Matrix(Matrix(:,2)==I,5)); %mean normal scores
        end
        s2=sum(Matrix(:,6))/(L-1);
        Z=realsqrt(sum(n.*Zbar.^2)/s2); %Statistics
        clear Matrix I
     case 'Wilcoxon test (paired samples)'
        x1=x(x(:,2)==1); x2=x(x(:,2)==2);
        dff=x2-x1; %difference between x1 and x2
        clear x1 x2
        nodiff = find(dff == 0);  %find null variations
        dff(nodiff) = []; %eliminate null variations
        if isempty(nodiff)==0 %tell me if there are null variations
            fprintf('There are %d null variations that will be deleted\n',length(nodiff))
            disp(' ')
        end
        if isempty(dff) %if all variations are null variations exit function
            disp('There are not variations. Van der Waerden - Wilcoxon test can''t be performed')
            return
        end
        clear nodiff %clear unnecessary variable
        %Ranks of absolute value of samples differences with sign
        R=tiedrank(abs(dff)); %ranks of diff
        P=(1+(R-c)./(length(R)+1-2*c))/2; %transform ranks in quantile
        Zp=norminv(P); %normal scores
        A=sign(dff).*Zp; %normal scores with sign
        Z=sum(A)/realsqrt(sum(A.^2)); %Statistics
        clear R P Zp A dff
    case 'Friedman test (two way ANOVA; without replication)'
        b=n(1); k=max(x(:,2)); %blocks and treatments
        y=reshape(x(:,1),b,k);
        Matrix=zeros(b,k); Rij=Matrix;
        for I=1:b
            Rij(I,:)=(tiedrank(y(I,:))); %ranks
        end
        Matrix=norminv((Rij-c)./(k+1-2*c)); %normal scores
        s2=sum(sum(Matrix.^2));
        Z=realsqrt(sum(sum(Matrix).^2)*(k-1)/s2);
        clear y Matrix I 
end

p=(1-normcdf(Z));  %p-value

%display results
tr=repmat('-',1,80); %set the divisor
fprintf('%s VERSION OF ',upper(method{type}))
switch avtest{tst}
    case 'Levene'
        fprintf('LEVENE TEST\n')
    case 'Mann-Whitney-Wilcoxon (unpaired samples)' 
        fprintf('MANN-WHITNEY-WILCOXON TEST\n')
    case 'Kruskal-Wallis (one way ANOVA)' 
        fprintf('KRUSKAL-WALLIS TEST\n')
    case 'Wilcoxon test (paired samples)'
        fprintf('WILCOXON TEST\n')
    case 'Friedman test (two way ANOVA; without replication)'
        fprintf('FRIEDMAN TEST FOR IDENTICAL TREATMENT EFFECTS:\n')
        disp('TWO-WAY BALANCED, COMPLETE BLOCK DESIGNS')
        disp(tr)
        disp(table(b*k,b,k,'VariableNames',{'Observations','Blocks','Treatments'}))
end
disp(tr)
disp('NORMAL SCORES STATISTICS')
disp(tr)
disp(table(Z,p,2*p,'VariableNames',{'Z','one_tailed_p_value','two_tailed_p_value'}))
disp(tr)

if nargout
    STATS.Z=Z;
    STATS.pvalue=[p 2*p];
end

ind=find(ismember({'Kruskal-Wallis (one way ANOVA)','Friedman test (two way ANOVA; without replication)'},avtest{tst}));

if ~isempty(ind) && 2*p<alpha
    disp(' ')
    disp('POST-HOC MULTIPLE COMPARISONS')
    disp(tr)
    switch ind
        case 1
            C=0.5*k*(k-1);
            alpha=alpha/C; %Bonferroni correction
            CD=tinv(1-alpha/2,L-k)*realsqrt(s2*(L-1-Z^2)/(L-k));
            Zdiff=zeros(k,k); mc=Zdiff; CDc=Zdiff;
            for J=1:k-1
                for I=J+1:k
                    Zdiff(I,J)=abs(Zbar(I)-Zbar(J));
                    CDc(I,J)=CD*realsqrt(1/sum(x(x(:,2)==J))+1/sum(x(x(:,2)==I)));
                    if Zdiff(I,J)>CDc(I,J)
                        mc(I,J)=1;
                    end
                end
            end
            %display results
            disp('Absolute difference among normal scores')
            disp(Zdiff)
            disp('Critical Values')
            disp(CDc)
            disp('Absolute difference > Critical Value')
            disp(mc)
        case 2
            tmp=repmat(sum(Rij),k,1); Rdiff=abs(tmp-tmp'); %Generate a matrix with the absolute differences among ranks
            gl=(b-1)*(k-1);
            cv=tinv(1-alpha/2,gl)*realsqrt((2*b*s2/gl)*(1-Z.^2/(b*(k-1)))); %critical value
            mc=Rdiff>cv; %Find differences greater than critical value
            %display results
            fprintf('Critical value: %0.4f\n',cv)
            disp('Absolute difference among mean ranks')
            disp(tril(Rdiff))
            disp('Absolute difference > Critical Value')
            disp(tril(mc))
    end
    if nargout
        STATS.posthoc=mc;
    end
end