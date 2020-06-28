#' fit a Bayesian model for the integration of nodes, interaction terms and networks
#'
#' One of the main functions in the JNNI package. Fits a path of JNNI models over different values of the tunning parameters.
#'
#' @param X a matrix of predictor variable
#' @param pathway a vector of pathway information, for instance, value 1 represent the corresponding node location belongs to pathway 1
#' @param L a list of Laplace information, each of which is a matrix capturing the relationship between node in the corresponding pathway, default is NULL
#' @param t a vector of the response variable
#' @param s_init a vector of tunning parameter related to the variance
#' @param a_init a vector of tunning parameter related to the initial selecting probability
#'
#' @return a list, each element corresponding to its tunning parameters contains values  as follows:
#' \item{mean_w_red}{The coefficients vector}
#' \item{beta_interaction}{The coefficients of interaction terms}
#' \item{beta_main}{The coefficients of main effects}
#' \item{BIC}{Bayesian Information Criterion}
#' \item{df}{The number of nonzero coefficients}
#' \item{index_path}{Selected networks' number}
#' @export
#'

#' @examples
#' set.seed(11)
#' u=runif(4,0.8,1.2)
#' u=round(u,2)
#' rrr=5
#' rate=round(1/sqrt(rrr),2)
#' n=300
#' p=1000
#' rho=0.6
#' a1=c(u[1],rate*rep(u[1],5))
#' a2=c(u[2],rate*rep(u[2],5))
#' a3=c(u[3],rate*rep(u[3],5))
#' b1=rate*rep(u[1],5)
#' b2=rate*rep(u[2],5)
#' b3=rate*rep(u[3],2)
#' c1=rate*rep(u[4],5)
#' M=100
#' set.seed(1)
#' X=matrix(NA,n,p)
#' p0=p/M
#' for(i in 1:n){
#' trans=rnorm(M)
#' for(j in 1:M){
#' start = (j-1)*p0+1
#' X[i,start]= trans[j]
#' xx=matrix(rnorm((p0-1), mean =(trans[j]*rho), sd=sqrt(1-rho^2)),ncol=1)
#' X[i,(start+1):(start+(p0-1))]=xx
#' }
#' }
#' x=list()
#' for(i in 1:M){
#' x[[i]]=X[,(1+(i-1)*p0):(p0*i)]
#' }
#'
#' Mm=M+M*(M-1)/2
#' Phix=list()
#' for(m in 1:M){
#' r=0
#' p_hat=p0+p0*(p0-1)/2
#' temp1=matrix(NA,n,p_hat)
#' temp1[,1:p0]=x[[m]]
#' for(i in 1:(p0-1)){
#' for(j in (1+i):p0){
#' r=r+1
#' temp1[,p0+r]=(x[[m]][,i])*(x[[m]][,j])
#' }
#' }
#' Phix[[m]]=temp1
#' }
#' cp=0
#' R1=matrix(0,M,M)
#' for(k in 1:(M-1)){
#' for(m in (k+1):M){
#' cp=cp+1
#' R1[k,m]=cp
#' temp1=matrix(NA,n,(p0*p0))
#' r=1
#' for(i in 1:p0){
#' temp1[,r:(r+p0-1)]=(x[[m]][,1:p0])*(x[[k]][,i])
#' r=r+p0
#' }
#' Phix[[M+cp]]=temp1
#' }
#' }
#'
#' real_w=list()
#' for(i in 1:M)
#' real_w[[i]]=rep(0,times=p0+(p0-1)*p0/2)
#' for(i in (1+M):Mm){
#' indexp1=which(R1==(i-M),arr.ind = TRUE)[1]
#' indexp2=which(R1==(i-M),arr.ind = TRUE)[2]
#' real_w[[i]]=rep(0,times=p0*p0)
#' }
#' real_w[[1]][1:length(a1)]=a1
#' real_w[[1]][(p0+1):(p0+length(b1))]=b1
#' real_w[[2]][1:length(a2)]=a2
#' real_w[[2]][(p0+1):(p0+length(b2))]=b2
#' real_w[[3]][1:length(a3)]=a3
#' real_w[[3]][(p0+1):(p0+length(b3))]=b3
#' real_w[[M+1]][1:length(c1)]=c1
#'
#'
#' real_y=0
#' for(m in c(1:3,(M+1))){
#' real_y=Phix[[m]]%*%real_w[[m]]+real_y
#' }
#' t = real_y + rnorm(n)
#' adj=matrix(0,p0,p0)
#' for(i in 2:p0)
#' adj[1,i]=1
#' adj = adj +t(adj)
#' C=NULL
#' for(k in 1:M)
#' C[[k]] = adj
#'
#' #laplacian matrix
#' L=NULL
#' for(k in 1:M){
#' dd = rowSums(adj)
#' ee = dd^(-0.5)
#' ee = diag(ee)
#' L[[k]] = diag(1,p0,p0)-ee%*%C[[k]]%*%ee
#' }
#'
#' pathway=matrix(0,p,1)
#' p_each=p0
#' for (i in 1:M){
#' pathway[(p_each*(i-1)+1):(p_each*i),1]=i
#' }
#' s_init=2e-4
#' a_init=0.96
#' ls=length(s_init)
#' la=length(a_init)
#' temp_1<-GANAIr.path(X,pathway,L=L,t,s_init,a_init)
GANAIr.path<-function(X,pathway,L=NULL,t,s_init,a_init){#######L is not including interacton
  n=dim(X)[1]
  n_var=dim(X)[2]
  M0=max(pathway)
  p_total=matrix(0,M0,1)
  for(i in 1:M0){
    p_total[i]=sum(pathway==i)
  }
  Tt=sum(t^2)

  p_int=p_total
  q_int=p_int*(p_int-1)/2
  M_int=M0
  p_sum_int=sum(p_int)
  p_hat_int=p_int*(p_int-1)/2+p_int
  pp_int=sum(p_hat_int[1:M0])
  b_s_int=p_sum_int+p_sum_int*(p_sum_int-1)/2


  if(is.null(L)){
    L=list()
    for(k in 1:M0){
      if(p_total[k]!=1){
        adj=matrix(1,p_total[k],p_total[k])
        adj=adj-diag(1,p_total[k],p_total[k])
        dd = rowSums(adj)
        ee = dd^(-0.5)
        ee = diag(ee)
        L[[k]] = diag(1,p_total[k],p_total[k])-ee%*%adj%*%ee
      }else{
        L[[k]]=1
      }
    }
  }


  LL_int=list()
  for(i in 1:M0){
    if((p_total[i]!=1)){
      LL_int[[i]]=matrix(0,p_hat_int[i],p_hat_int[i])
      LL_int[[i]][1:p_total[i],1:p_total[i]]=L[[i]]
      LL_int[[i]][(p_total[i]+1):p_hat_int[i],(p_total[i]+1):p_hat_int[i]]=temp_function2(L[[i]])
    }else{
      LL_int[[i]]=L[[i]]
    }
  }

  L1_int=matrix(0,pp_int,pp_int)
  s=1
  for (m in 1:M_int){
    L1_int[s:(s+p_hat_int[m]-1),s:(s+p_hat_int[m]-1)]=LL_int[[m]]
    s=s+p_hat_int[m]
  }

  R_int=lapply(p_int,location)


  node_location_int=matrix(1,p_sum_int,1)
  node_location_int[1:p_int[1]]=1:p_int[1]
  for (m in 2:M_int){
    node_location_int[(sum(p_int[1:(m-1)])+1):(sum(p_int[1:(m-1)])+p_int[m])]<-(sum(p_hat_int[1:(m-1)])+1):(sum(p_hat_int[1:(m-1)])+p_int[m])
  }
  within_interaction_location_int=setdiff(1:pp_int,node_location_int)
  within_interaction_location2_int=matrix(FALSE,p_sum_int,p_sum_int)
  within_interaction_location2_int[1:p_int[1],1:p_int[1]]=TRUE
  for (m in 2:M_int){
    within_interaction_location2_int[(sum(p_int[1:(m-1)])+1):(sum(p_int[1:(m)])),(sum(p_int[1:(m-1)])+1):(sum(p_int[1:(m)]))]=TRUE
  }


  out_interaction_location2_int=!(within_interaction_location2_int)
  RR_int=matrix(0,p_sum_int,p_sum_int)
  temp=(lower.tri(out_interaction_location2_int))*out_interaction_location2_int
  RR_int[which(temp==1)]=c(1:(b_s_int-pp_int))

  Phi_int=matrix(0,n,pp_int)
  Phi_int[,node_location_int]=X

  s=0
  r1=1
  for(m in 1:M_int){
    if(q_int[m]>0){
      temp=matrix(0,n,q_int[m])
      r=0
      for(i in 1:(p_int[m]-1)){
        temp1=p_int[m]-i
        temp[,(r+1):(r+temp1)]=(X[,(s+i)])*(X[,(s+i+1):(p_int[m]+s)])
        r=r+temp1
      }
      Phi_int[,(r1+p_int[m]):(r1+p_hat_int[m]-1)]=temp
      s=s+p_int[m]
      r1=r1+p_hat_int[m]
    }
  }


  r=1
  temp2=0
  s=0
  Phi2_int=matrix(0,n,(b_s_int-pp_int))
  for(m in 1:(M_int-1)){
    temp=p_sum_int-temp2-p_int[m]
    for(i in 1:p_int[m]){
      Phi2_int[,((i-1)*temp+1+s):(i*temp+s)]=X[,(temp2+i)]*X[,(p_int[m]+temp2+1):p_sum_int]
    }
    temp2=temp2+p_int[m]
    s=temp*p_int[m]+s
  }



  Xt_int=matrix(0,b_s_int,1)
  Xt_int[1:pp_int,]=t(Phi_int)%*%t
  Xt_int[(pp_int+1):b_s_int,]=t(Phi2_int)%*%t
  cols_PP1_int=colSums(Phi_int*Phi_int)
  cols_PP2_int=colSums((Phi2_int*Phi2_int))
  pathway_new_int=matrix(0,pp_int,1)
  pathway_new_int[node_location_int]=pathway


  ls=length(s_init)
  la=length(a_init)
  result_temp=list()
  lll=1
  s1=1
  for(i in 1:ls){
    s2=s_init[i]
    print(paste0('The ',i,'th s'))
    for(j in 1:la){
      a=a_init[j]
      p=p_int
      M=M_int
      p_sum=p_sum_int
      p_hat=p_hat_int
      pp=pp_int
      b_s=b_s_int
      L1=L1_int
      R=R_int
      node_location=node_location_int
      within_interaction_location=within_interaction_location_int
      within_interaction_location2=within_interaction_location2_int
      out_interaction_location2=out_interaction_location2_int
      RR=RR_int
      Phi=Phi_int
      q=q_int
      Phi2=Phi2_int
      Xt=Xt_int
      cols_PP1=cols_PP1_int
      cols_PP2= cols_PP2_int
      gamma=a*rep(1,times=M)
      eta_ori=matrix(a,b_s,1)
      b_hat= rep(1,times=M)
      a_hat = rep(1,times=M)
      # for tau
      g0 = exp(-12)
      h0 = exp(-12)
      a=1
      b=1
      c=1
      d=1
      mean_w_ori=rep(0,times=b_s)
      d_hat_ori=c_hat_ori=rep(1,times=b_s)
      y1=t
      tau =as.numeric(sqrt(1/var(t)))
      niter = 0
      tol=1e-3
      diff_sum=1
      index_path_ori=c(1:M)

      sigma_vec=sigma_vec_ori=matrix(0,b_s,1)
      mean_w=mean_w_ori
      pathway_inter_includ=matrix(0,pp,1)
      s=1
      for (i in 1:M){
        pathway_inter_includ[s:(s+p_hat[i]-1),1]=i
        s=s+p_hat[i]
      }
      locat_within=c(1:pp)
      eta=eta_ori

      while (diff_sum>tol&niter<1000){
        niter = niter+1

        rm(Phi2)
        gamma_new=matrix(0,pp,pp)
        gamma_new[1:p_hat[1],1:p_hat[1]]=gamma[1]
        if(M>1){
          for (m in 2:M){
            gamma_new[(sum(p_hat[1:(m-1)])+1):(sum(p_hat[1:m])),(sum(p_hat[1:(m-1)])+1):(sum(p_hat[1:m]))]=gamma[m]
          }
        }
        eta_temp=matrix(eta[node_location],p_sum,1)
        temp=eta_temp%*%t(eta_temp)
        temp2=as.vector(t(upper.tri(temp)*within_interaction_location2))
        eta_prod=(as.vector(temp))[temp2==1]


        A=gamma_new/s1*L1+(1-gamma_new)/s2*diag(1,pp,pp)
        rm(gamma_new)
        B=matrix(0,pp,pp)
        B[node_location,node_location]=diag(eta[node_location],length(node_location),length(node_location))/s1+diag(1-eta[node_location],length(node_location),length(node_location))/s2
        B[within_interaction_location,within_interaction_location]=diag((eta[within_interaction_location]+eta_prod)/s1+(2-eta[within_interaction_location]-eta_prod)/s2,length(within_interaction_location),length(within_interaction_location))
        A=A+B
        rm(B)

        sigma_vec[1:pp]=1/(tau*cols_PP1[1:pp]+(diag(A))[1:pp])


        mean_w_old=mean_w_ori

        for(i in 1:pp){
          sum_a=sum(A[,i]*mean_w[1:pp])-A[i,i]*mean_w[i]
          y1=y1+Phi[,i]*mean_w[i]
          mean_w[i]=sigma_vec[i]*(sum(Phi[,i]*y1)*tau-0.5*sum_a)
          y1=y1-Phi[,i]*mean_w[i]
        }
        wwmtp_path_within=diag(sigma_vec[1:pp],pp,pp)+mean_w[1:pp]%*%t(mean_w[1:pp])
        sigma_vec_ori[locat_within]=sigma_vec[1:pp]
        rm(A)

        mean_w_ori[1:pp_int]=rep(0,pp_int)
        mean_w_ori[locat_within]=mean_w[1:pp]


        #update Q(gamma)
        s=1
        gamma=rep(0,M)
        for(m in 1:M){
          temp=digamma(b_hat[m])-digamma(a_hat[m])+0.5*p_hat[m]*log(s1/s2)+-0.5*log(det(L1[(s:(s+p_hat[m]-1)),(s:(s+p_hat[m]-1))]+1e-6*diag(1,p_hat[m],p_hat[m])))+0.5*sum(diag(wwmtp_path_within[s:(s+p_hat[m]-1),s:(s+p_hat[m]-1)]%*%(L1[s:(s+p_hat[m]-1),s:(s+p_hat[m]-1)]/s1-diag(1,p_hat[m])/s2)))
          gamma[m]= 1/(1+exp(temp))
          s=p_hat[m]+s
        }
        rm(wwmtp_path_within)

        #selected pathway interaction updating
        gamma_trunc=ifelse(gamma>0.5,1,0)
        index_path1=which(gamma>0.5)
        if(length(index_path1)==0) break;
        index_path_ori=gamma_trunc*(index_path_ori[index_path_ori!=0])#######

        #update Q(u)
        a_hat = (a+ gamma)[index_path1]
        b_hat = (b + 1 - gamma)[index_path1]

        p=p_total[index_path_ori]
        M=length(p)
        p_sum=sum(p)
        q=p*(p-1)/2
        p_hat=p*(p-1)/2+p
        pp=sum(p_hat[1:M])
        b_s=p_sum+p_sum*(p_sum-1)/2
        pathway_second=matrix(0,p_sum,1)
        s=1
        for (i in 1:M){
          pathway_second[s:(s+p[i]-1),1]=i
          s=s+p[i]
        }

        locat_within=which(is.element(pathway_inter_includ,index_path_ori))
        temp=which(is.element(pathway,index_path_ori))
        temp1=RR_int[temp,temp]
        locat_inter_out=as.vector(temp1[temp1!=0])+pp_int


        mean_w=mean_w_ori[c(locat_within,locat_inter_out)]
        eta=eta_ori[c(locat_within,locat_inter_out)]
        c_hat=c_hat_ori[c(locat_within,locat_inter_out)]
        d_hat=d_hat_ori[c(locat_within,locat_inter_out)]
        LL=LL_int[index_path_ori]

        L1=matrix(0,pp,pp)
        s=1
        for (m in 1:M){
          L1[s:(s+p_hat[m]-1),s:(s+p_hat[m]-1)]=LL[[m]]
          s=s+p_hat[m]
        }

        R=lapply(p,location)

        node_location=matrix(1,p_sum,1)
        node_location[1:p[1]]=1:p[1]
        if(M>1){
          for (m in 2:M){
            node_location[(sum(p[1:(m-1)])+1):(sum(p[1:(m-1)])+p[m])]<-(sum(p_hat[1:(m-1)])+1):(sum(p_hat[1:(m-1)])+p[m])
          }
        }
        within_interaction_location=setdiff(1:pp,node_location)
        within_interaction_location2=matrix(FALSE,p_sum,p_sum)
        within_interaction_location2[1:p[1],1:p[1]]=TRUE
        if(M>1){
          for (m in 2:M){
            within_interaction_location2[(sum(p[1:(m-1)])+1):(sum(p[1:(m)])),(sum(p[1:(m-1)])+1):(sum(p[1:(m)]))]=TRUE
          }
        }
        out_interaction_location2=!(within_interaction_location2)
        RR=matrix(0,p_sum,p_sum)
        temp=(lower.tri(out_interaction_location2))*out_interaction_location2
        RR[which(temp==1)]=c(1:(b_s-pp))


        Phi=matrix(Phi_int[,locat_within],n,pp)

        Xt=matrix(0,b_s,1)
        Xt[1:pp,]=Xt_int[locat_within]
        Phi2=matrix(0,n,(b_s-pp))
        if(M>1){
          Phi2=as.matrix(Phi2_int[,locat_inter_out-pp_int])
          Xt[(pp+1):b_s,]=Xt_int[locat_inter_out]
        }




        cols_PP1=colSums(Phi*Phi)
        cols_PP2=colSums((Phi2*Phi2))
        pathway_new=matrix(0,pp,1)
        pathway_new[node_location]=pathway_second

        sigma_vec=rep(0,b_s)

        #interaction updating
        eta_temp=matrix(eta[node_location],p_sum,1)
        temp=eta_temp%*%t(eta_temp)
        temp2=as.vector(t(upper.tri(temp)*within_interaction_location2))
        eta_prod=(as.vector(temp))[temp2==1]
        temp2=as.vector(t(upper.tri(temp)*out_interaction_location2))
        eta_prod_out=(as.vector(temp))[temp2==1]
        summ=(eta[(pp+1):b_s]+eta_prod_out)/s1+(2-eta[(pp+1):b_s]-eta_prod_out)/s2
        if(M>1){
          sigma_vec[(pp+1):b_s]=as.vector(1/(tau*cols_PP2+summ))
          y1=t-Phi%*%mean_w[1:pp]-Phi2%*%mean_w[(pp+1):b_s]
          for(i in (1+pp):b_s){
            y1=y1+Phi2[,(i-pp)]*mean_w[i]
            mean_w[i]=(sum(Phi2[,(i-pp)]*y1))*tau*sigma_vec[i]
            y1=y1-Phi2[,(i-pp)]*mean_w[i]
          }
        }else{
          y1=t-Phi%*%mean_w[1:pp]
        }


        mean_w_ori[(pp_int+1):b_s_int]=rep(0,(b_s_int-pp_int))
        if(M>1){
          mean_w_ori[locat_inter_out]=mean_w[(pp+1):b_s]
          wwmtp_path_out=sigma_vec[(pp+1):b_s]+mean_w[(pp+1):b_s]^2
        }else{
          wwmtp_path_out=0
        }


        wwmtp_path_within=diag(sigma_vec_ori[locat_within],length(locat_within),length(locat_within))+mean_w[1:pp]%*%t(mean_w[1:pp])

        for(i in node_location){
          m=pathway_new[i]
          id_node=which(node_location==i)
          if((p[m]!=1)){
            id_temp=as.vector(setdiff(which(pathway_new==pathway_new[i]),i))
            id_nodep=ifelse(m>1,i-sum(p_hat[1:(m-1)]),i)
            id_temp1=ifelse(rep(m,(p[m]-1))>1,id_temp-sum(p_hat[1:(m-1)]),id_temp)
            temp=pmax((R[[m]][id_nodep,id_temp1]),(R[[m]][id_temp1,id_nodep]))
            temp=ifelse(rep(m,(p[m]-1))>1,temp+sum(p_hat[1:(m-1)]),temp)+p[m]
            sum_kj=sum(eta[id_temp]*(0.5*log(s1/s2)+0.5*diag(wwmtp_path_within[temp,temp])*(1/s1-1/s2)))
          }else{
            sum_kj=0
          }
          id_temp=as.vector(setdiff(1:p_sum,which(pathway_second==m)))
          id_temp1=as.vector(which((pathway_new!=m)&(pathway_new!=0)))
          temp=pmax(RR[id_node,id_temp],RR[id_temp,id_node])
          sum_mkj=sum(eta[id_temp1]*(0.5*log(s1/s2)+0.5*wwmtp_path_out[temp]*(1/s1-1/s2)))
          eta[i]=1/(1 +exp(0.5*log(s1/s2)+ 0.5*wwmtp_path_within[i,i]*(1/s1-1/s2)+ digamma(d_hat[i])-digamma(c_hat[i])+sum_kj+sum_mkj))
        }
        eta[within_interaction_location]=1/(1 +exp(0.5*log(s1/s2)+ 0.5*(diag(wwmtp_path_within)[within_interaction_location])*(1/s1-1/s2)+ digamma(d_hat[within_interaction_location])-digamma(c_hat[within_interaction_location])))

        if(M>1) eta[(pp+1):b_s]=1/(1 +exp(0.5*log(s1/s2)+ 0.5*wwmtp_path_out[1:(b_s-pp)]*(1/s1-1/s2)+digamma(d_hat[(pp+1):b_s])-digamma(c_hat[(pp+1):b_s])))
        rm(wwmtp_path_within)

        #update Q(v) for each pathway
        c_hat= c + eta
        d_hat= d + 1 - eta


        w_mul_Phi=t-y1
        trace0=sum(w_mul_Phi^2)
        trace0=trace0+sum(cols_PP1*sigma_vec_ori[locat_within])
        if(M>1) trace0=trace0+sum(cols_PP2*sigma_vec[(pp+1):b_s])
        mean_w_Xt=sum(mean_w*Xt)


        g_hat = g0 + n/2
        h_hat = h0 + 0.5*Tt - mean_w_Xt + 0.5*trace0

        #update Q(tau)
        tau = g_hat/h_hat

        eta_ori=rep(0,b_s_int)
        eta_ori[locat_within]=eta[1:pp]
        if(M>1){
          eta_ori[locat_inter_out]=eta[(pp+1):b_s]
          c_hat_ori[locat_inter_out]=c_hat[(pp+1):b_s]
          d_hat_ori[locat_inter_out]=d_hat[(pp+1):b_s]
        }
        c_hat_ori[locat_within]=c_hat[1:pp]
        d_hat_ori[locat_within]=d_hat[1:pp]


        temp=which((mean_w_old-mean_w_ori)!=0)

        diff_sum = mean(abs((mean_w_old-mean_w_ori)[temp]/mean_w_old[temp]))
      }

      if(length(index_path1)==0) {
        eta_ori=mean_w_ori=rep(0,b_s_int)
        eta=mean_w=rep(0,length(eta))
      }
      index_pred1=which(eta>0.5)
      temp=rep(0,b_s)
      temp[index_pred1]=mean_w[index_pred1]
      if(M>1){
        error=mean((t-Phi%*%temp[1:pp]-Phi2%*%temp[(pp+1):b_s])^2)
      }else{
        error=mean((t-Phi%*%temp[1:pp])^2)
      }




      index_pred=which(eta_ori>0.5)
      mean_w_red=rep(0,b_s_int)
      mean_w_red[index_pred]=mean_w_ori[index_pred]
      df=length(which(eta>0.5))
      BIC=log(error)+df*log(n)/n
      beta_interaction=mean_w_red[c(within_interaction_location_int,c((pp_int+1):b_s_int))]
      beta_main=mean_w_red[node_location_int]
      temp_2=list(mean_w_red=mean_w_red,beta_interaction=beta_interaction,beta_main=beta_main,BIC=BIC,df=df,index_path=index_path_ori)
      result_temp[[lll]]<-temp_2

      ##just for observe
      cat(paste0('the',lll,'th:'),'df=',df,'\n')
      lll=lll+1

    }
  }

  return(result_temp)

}




#' fit a Bayesian model for the integration of nodes, interaction terms and networks
#'
#' One of the main functions in the JNNI package.
#'
#' @param X a matrix of predictor variable
#' @param pathway a vector of pathway information, for instance, value 1 represent the corresponding node location belongs to pathway 1
#' @param L a list of Laplace information, each of which is a matrix capturing the relationship between node in the corresponding pathway, default is NULL
#' @param t a vector of the response variable
#' @param s2 a numeric
#' @param prob a numeric
#'
#' @return A list containing values as follows:
#' \item{mean_w_red}{The coefficients vector}
#' \item{beta_interaction}{The coefficients of interaction terms}
#' \item{beta_main}{The coefficients of main effects}
#' \item{BIC}{Bayesian Information Criterion}
#' \item{df}{The number of nonzero coefficients}
#' \item{index_path}{Selected networks' number}
#' @export
#'
#' @examples
#' library(JNNI)
#' set.seed(1)
#' X=matrix(rnorm(50),10,5)
#' pathway=c(1,1,2,2,3)
#' b=c(3,3,0,0,0)
#' t=as.numeric(X%*%b+rnorm(10))
#' s2=1e-3
#' prob=0.96
#' temp=GANAIr(X,pathway,L=NULL,t,s2,prob)
GANAIr<-function(X,pathway,L=NULL,t,s2,prob){
  n=dim(X)[1]
  n_var=dim(X)[2]
  M0=max(pathway)
  p_total=matrix(0,M0,1)
  for(i in 1:M0){
    p_total[i]=sum(pathway==i)
  }
  Tt=sum(t^2)

  p=p_total
  q=p*(p-1)/2
  M=M0
  p_sum=sum(p)
  p_hat=p*(p-1)/2+p
  pp_int=sum(p_hat[1:M0])
  b_s_int=p_sum+p_sum*(p_sum-1)/2

  if(is.null(L)){
    L=list()
    for(k in 1:M0){
      if(p_total[k]!=1){
        adj=matrix(1,p_total[k],p_total[k])
        adj=adj-diag(1,p_total[k],p_total[k])
        dd = rowSums(adj)
        ee = dd^(-0.5)
        ee = diag(ee)
        L[[k]] = diag(1,p_total[k],p_total[k])-ee%*%adj%*%ee
      }else{
        L[[k]]=1
      }
    }
  }

  LL_int=list()
  for(i in 1:M0){
    if((p_total[i]!=1)){
      LL_int[[i]]=matrix(0,p_hat[i],p_hat[i])
      LL_int[[i]][1:p_total[i],1:p_total[i]]=L[[i]]
      LL_int[[i]][(p_total[i]+1):p_hat[i],(p_total[i]+1):p_hat[i]]=temp_function2(L[[i]])##diag(1,q_int[i],q_int[i])
    }else{
      LL_int[[i]]=L[[i]]
    }
  }

  L1=matrix(0,pp_int,pp_int)
  s=1
  for (m in 1:M){
    L1[s:(s+p_hat[m]-1),s:(s+p_hat[m]-1)]=LL_int[[m]]
    s=s+p_hat[m]
  }

  R=lapply(p,location)

  node_location_int=matrix(1,p_sum,1)
  node_location_int[1:p[1]]=1:p[1]
  for (m in 2:M){
    node_location_int[(sum(p[1:(m-1)])+1):(sum(p[1:(m-1)])+p[m])]<-(sum(p_hat[1:(m-1)])+1):(sum(p_hat[1:(m-1)])+p[m])
  }
  within_interaction_location_int=setdiff(1:pp_int,node_location_int)
  within_interaction_location2=matrix(FALSE,p_sum,p_sum)
  within_interaction_location2[1:p[1],1:p[1]]=TRUE
  for (m in 2:M){
    within_interaction_location2[(sum(p[1:(m-1)])+1):(sum(p[1:(m)])),(sum(p[1:(m-1)])+1):(sum(p[1:(m)]))]=TRUE
  }


  out_interaction_location2=!(within_interaction_location2)
  RR_int=matrix(0,p_sum,p_sum)
  temp=(lower.tri(out_interaction_location2))*out_interaction_location2
  RR_int[which(temp==1)]=c(1:(b_s_int-pp_int))

  Phi_int=matrix(0,n,pp_int)
  Phi_int[,node_location_int]=X

  s=0
  r1=1
  for(m in 1:M){
    if(q[m]>0){
      temp=matrix(0,n,q[m])
      r=0
      for(i in 1:(p[m]-1)){
        temp1=p[m]-i
        temp[,(r+1):(r+temp1)]=(X[,(s+i)])*(X[,(s+i+1):(p[m]+s)])
        r=r+temp1
      }
      Phi_int[,(r1+p[m]):(r1+p_hat[m]-1)]=temp
      s=s+p[m]
      r1=r1+p_hat[m]
    }
  }


  r=1
  temp2=0
  s=0
  Phi2_int=matrix(0,n,(b_s_int-pp_int))
  for(m in 1:(M-1)){
    temp=p_sum-temp2-p[m]
    for(i in 1:p[m]){
      Phi2_int[,((i-1)*temp+1+s):(i*temp+s)]=X[,(temp2+i)]*X[,(p[m]+temp2+1):p_sum]
    }
    temp2=temp2+p[m]
    s=temp*p[m]+s
  }



  Xt_int=matrix(0,b_s_int,1)
  Xt_int[1:pp_int,]=t(Phi_int)%*%t
  Xt_int[(pp_int+1):b_s_int,]=t(Phi2_int)%*%t


  cols_PP1=colSums(Phi_int*Phi_int)
  cols_PP2=colSums((Phi2_int*Phi2_int))
  pathway_new=matrix(0,pp_int,1)
  pathway_new[node_location_int]=pathway


  result_temp=list()
  s1=1


  a=prob
  pp=pp_int
  b_s=b_s_int
  node_location=node_location_int
  within_interaction_location=within_interaction_location_int
  RR=RR_int
  Phi=Phi_int
  Phi2=Phi2_int
  Xt=Xt_int
  gamma=a*rep(1,times=M)
  eta_ori=matrix(a,b_s,1)
  b_hat= rep(1,times=M)
  a_hat = rep(1,times=M)
  # for tau
  g0 = exp(-12)
  h0 = exp(-12)
  a=1
  b=1
  c=1
  d=1
  mean_w_ori=rep(0,times=b_s)
  d_hat_ori=c_hat_ori=rep(1,times=b_s)
  y1=t
  tau =as.numeric(sqrt(1/var(t)))
  niter = 0
  tol=1e-3
  diff_sum=1
  index_path_ori=c(1:M)

  sigma_vec=sigma_vec_ori=matrix(0,b_s,1)
  mean_w=mean_w_ori
  pathway_inter_includ=matrix(0,pp,1)
  s=1
  for (i in 1:M){
    pathway_inter_includ[s:(s+p_hat[i]-1),1]=i
    s=s+p_hat[i]
  }
  locat_within=c(1:pp)
  eta=eta_ori

  while (diff_sum>tol&niter<1000){
    niter = niter+1
    gamma_new=matrix(0,pp,pp)
    gamma_new[1:p_hat[1],1:p_hat[1]]=gamma[1]
    if(M>1){
      for (m in 2:M){
        gamma_new[(sum(p_hat[1:(m-1)])+1):(sum(p_hat[1:m])),(sum(p_hat[1:(m-1)])+1):(sum(p_hat[1:m]))]=gamma[m]
      }
    }
    eta_temp=matrix(eta[node_location],p_sum,1)
    temp=eta_temp%*%t(eta_temp)
    temp2=as.vector(t(upper.tri(temp)*within_interaction_location2))
    eta_prod=(as.vector(temp))[temp2==1]
    A=gamma_new/s1*L1+(1-gamma_new)/s2*diag(1,pp,pp)
    rm(gamma_new)
    B=matrix(0,pp,pp)
    B[node_location,node_location]=diag(eta[node_location],length(node_location),length(node_location))/s1+diag(1-eta[node_location],length(node_location),length(node_location))/s2
    B[within_interaction_location,within_interaction_location]=diag((eta[within_interaction_location]+eta_prod)/s1+(2-eta[within_interaction_location]-eta_prod)/s2,length(within_interaction_location),length(within_interaction_location))
    A=A+B
    rm(B)

    sigma_vec[1:pp]=1/(tau*cols_PP1[1:pp]+(diag(A))[1:pp])
    mean_w_old=mean_w_ori

    for(i in 1:pp){
      sum_a=sum(A[,i]*mean_w[1:pp])-A[i,i]*mean_w[i]
      y1=y1+Phi[,i]*mean_w[i]
      mean_w[i]=sigma_vec[i]*(sum(Phi[,i]*y1)*tau-0.5*sum_a)
      y1=y1-Phi[,i]*mean_w[i]
    }
    wwmtp_path_within=diag(sigma_vec[1:pp],pp,pp)+mean_w[1:pp]%*%t(mean_w[1:pp])
    sigma_vec_ori[locat_within]=sigma_vec[1:pp]
    rm(A)

    mean_w_ori[1:pp_int]=rep(0,pp_int)
    mean_w_ori[locat_within]=mean_w[1:pp]


    #update Q(gamma)
    s=1
    gamma=rep(0,M)
    for(m in 1:M){
      temp=digamma(b_hat[m])-digamma(a_hat[m])+0.5*p_hat[m]*log(s1/s2)+-0.5*log(det(L1[(s:(s+p_hat[m]-1)),(s:(s+p_hat[m]-1))]+1e-6*diag(1,p_hat[m],p_hat[m])))+0.5*sum(diag(wwmtp_path_within[s:(s+p_hat[m]-1),s:(s+p_hat[m]-1)]%*%(L1[s:(s+p_hat[m]-1),s:(s+p_hat[m]-1)]/s1-diag(1,p_hat[m])/s2)))
      gamma[m]= 1/(1+exp(temp))
      s=p_hat[m]+s
    }

    gamma_trunc=ifelse(gamma>0.5,1,0)
    index_path1=which(gamma>0.5)
    if(length(index_path1)==0) break;
    index_path_ori=gamma_trunc*(index_path_ori[index_path_ori!=0])#######

    #update Q(u)
    a_hat = (a+ gamma)[index_path1]
    b_hat = (b + 1 - gamma)[index_path1]

    p=p_total[index_path_ori]
    M=length(p)
    p_sum=sum(p)
    q=p*(p-1)/2
    p_hat=p*(p-1)/2+p
    pp=sum(p_hat[1:M])
    b_s=p_sum+p_sum*(p_sum-1)/2
    pathway_second=matrix(0,p_sum,1)
    s=1
    for (i in 1:M){
      pathway_second[s:(s+p[i]-1),1]=i
      s=s+p[i]
    }

    locat_within=which(is.element(pathway_inter_includ,index_path_ori))
    temp=which(is.element(pathway,index_path_ori))
    temp1=RR_int[temp,temp]
    locat_inter_out=as.vector(temp1[temp1!=0])+pp_int


    mean_w=mean_w_ori[c(locat_within,locat_inter_out)]
    eta=eta_ori[c(locat_within,locat_inter_out)]
    c_hat=c_hat_ori[c(locat_within,locat_inter_out)]
    d_hat=d_hat_ori[c(locat_within,locat_inter_out)]
    LL=LL_int[index_path_ori]

    L1=matrix(0,pp,pp)
    s=1
    for (m in 1:M){
      L1[s:(s+p_hat[m]-1),s:(s+p_hat[m]-1)]=LL[[m]]
      s=s+p_hat[m]
    }

    R=lapply(p,location)

    node_location=matrix(1,p_sum,1)
    node_location[1:p[1]]=1:p[1]
    if(M>1){
      for (m in 2:M){
        node_location[(sum(p[1:(m-1)])+1):(sum(p[1:(m-1)])+p[m])]<-(sum(p_hat[1:(m-1)])+1):(sum(p_hat[1:(m-1)])+p[m])
      }
    }
    within_interaction_location=setdiff(1:pp,node_location)
    within_interaction_location2=matrix(FALSE,p_sum,p_sum)
    within_interaction_location2[1:p[1],1:p[1]]=TRUE
    if(M>1){
      for (m in 2:M){
        within_interaction_location2[(sum(p[1:(m-1)])+1):(sum(p[1:(m)])),(sum(p[1:(m-1)])+1):(sum(p[1:(m)]))]=TRUE
      }
    }
    out_interaction_location2=!(within_interaction_location2)
    RR=matrix(0,p_sum,p_sum)
    temp=(lower.tri(out_interaction_location2))*out_interaction_location2
    RR[which(temp==1)]=c(1:(b_s-pp))


    Phi=matrix(Phi_int[,locat_within],n,pp)

    Xt=matrix(0,b_s,1)
    Xt[1:pp,]=Xt_int[locat_within]
    Phi2=matrix(0,n,(b_s-pp))
    if(M>1){
      Phi2=as.matrix(Phi2_int[,locat_inter_out-pp_int])
      Xt[(pp+1):b_s,]=Xt_int[locat_inter_out]
    }




    cols_PP1=colSums(Phi*Phi)
    cols_PP2=colSums((Phi2*Phi2))


    pathway_new=matrix(0,pp,1)
    pathway_new[node_location]=pathway_second

    sigma_vec=rep(0,b_s)

    #interaction updating
    eta_temp=matrix(eta[node_location],p_sum,1)
    temp=eta_temp%*%t(eta_temp)
    temp2=as.vector(t(upper.tri(temp)*within_interaction_location2))
    eta_prod=(as.vector(temp))[temp2==1]
    temp2=as.vector(t(upper.tri(temp)*out_interaction_location2))
    eta_prod_out=(as.vector(temp))[temp2==1]
    summ=(eta[(pp+1):b_s]+eta_prod_out)/s1+(2-eta[(pp+1):b_s]-eta_prod_out)/s2
    if(M>1){
      sigma_vec[(pp+1):b_s]=as.vector(1/(tau*cols_PP2+summ))
      y1=t-Phi%*%mean_w[1:pp]-Phi2%*%mean_w[(pp+1):b_s]
      for(i in (1+pp):b_s){
        y1=y1+Phi2[,(i-pp)]*mean_w[i]
        mean_w[i]=(sum(Phi2[,(i-pp)]*y1))*tau*sigma_vec[i]
        y1=y1-Phi2[,(i-pp)]*mean_w[i]
      }
    }else{
      y1=t-Phi%*%mean_w[1:pp]
    }


    mean_w_ori[(pp_int+1):b_s_int]=rep(0,(b_s_int-pp_int))
    if(M>1){
      mean_w_ori[locat_inter_out]=mean_w[(pp+1):b_s]
      wwmtp_path_out=sigma_vec[(pp+1):b_s]+mean_w[(pp+1):b_s]^2
    }else{
      wwmtp_path_out=0
    }


    wwmtp_path_within=diag(sigma_vec_ori[locat_within],length(locat_within),length(locat_within))+mean_w[1:pp]%*%t(mean_w[1:pp])

    for(i in node_location){
      m=pathway_new[i]
      id_node=which(node_location==i)
      if((p[m]!=1)){
        id_temp=as.vector(setdiff(which(pathway_new==pathway_new[i]),i))
        id_nodep=ifelse(m>1,i-sum(p_hat[1:(m-1)]),i)

        id_temp1=ifelse(rep(m,(p[m]-1))>1,id_temp-sum(p_hat[1:(m-1)]),id_temp)
        temp=pmax((R[[m]][id_nodep,id_temp1]),(R[[m]][id_temp1,id_nodep]))
        temp=ifelse(rep(m,(p[m]-1))>1,temp+sum(p_hat[1:(m-1)]),temp)+p[m]
        sum_kj=sum(eta[id_temp]*(0.5*log(s1/s2)+0.5*diag(wwmtp_path_within[temp,temp])*(1/s1-1/s2)))
      }else{
        sum_kj=0
      }
      id_temp=as.vector(setdiff(1:p_sum,which(pathway_second==m)))
      id_temp1=as.vector(which((pathway_new!=m)&(pathway_new!=0)))
      temp=pmax(RR[id_node,id_temp],RR[id_temp,id_node])
      sum_mkj=sum(eta[id_temp1]*(0.5*log(s1/s2)+0.5*wwmtp_path_out[temp]*(1/s1-1/s2)))
      eta[i]=1/(1 +exp(0.5*log(s1/s2)+ 0.5*wwmtp_path_within[i,i]*(1/s1-1/s2)+ digamma(d_hat[i])-digamma(c_hat[i])+sum_kj+sum_mkj))
    }
    eta[within_interaction_location]=1/(1 +exp(0.5*log(s1/s2)+ 0.5*(diag(wwmtp_path_within)[within_interaction_location])*(1/s1-1/s2)+ digamma(d_hat[within_interaction_location])-digamma(c_hat[within_interaction_location])))

    if(M>1) eta[(pp+1):b_s]=1/(1 +exp(0.5*log(s1/s2)+ 0.5*wwmtp_path_out[1:(b_s-pp)]*(1/s1-1/s2)+digamma(d_hat[(pp+1):b_s])-digamma(c_hat[(pp+1):b_s])))
    rm(wwmtp_path_within)

    #update Q(v) for each pathway
    c_hat= c + eta
    d_hat= d + 1 - eta


    w_mul_Phi=t-y1
    trace0=sum(w_mul_Phi^2)
    trace0=trace0+sum(cols_PP1*sigma_vec_ori[locat_within])
    if(M>1) trace0=trace0+sum(cols_PP2*sigma_vec[(pp+1):b_s])
    mean_w_Xt=sum(mean_w*Xt)


    g_hat = g0 + n/2
    h_hat = h0 + 0.5*Tt - mean_w_Xt + 0.5*trace0

    #update Q(tau)
    tau = g_hat/h_hat

    eta_ori=rep(0,b_s_int)
    eta_ori[locat_within]=eta[1:pp]
    if(M>1){
      eta_ori[locat_inter_out]=eta[(pp+1):b_s]
      c_hat_ori[locat_inter_out]=c_hat[(pp+1):b_s]
      d_hat_ori[locat_inter_out]=d_hat[(pp+1):b_s]
    }
    c_hat_ori[locat_within]=c_hat[1:pp]
    d_hat_ori[locat_within]=d_hat[1:pp]


    temp=which((mean_w_old-mean_w_ori)!=0)

    diff_sum = mean(abs((mean_w_old-mean_w_ori)[temp]/mean_w_old[temp]))
  }

  if(length(index_path1)==0) {
    eta_ori=mean_w_ori=rep(0,b_s_int)
    eta=mean_w=rep(0,length(eta))
  }
  index_pred1=which(eta>0.5)
  temp=rep(0,b_s)
  temp[index_pred1]=mean_w[index_pred1]
  if(M>1){
    error=mean((t-Phi%*%temp[1:pp]-Phi2%*%temp[(pp+1):b_s])^2)
  }else{
    error=mean((t-Phi%*%temp[1:pp])^2)
  }




  index_pred=which(eta_ori>0.5)
  mean_w_red=rep(0,b_s_int)
  mean_w_red[index_pred]=mean_w_ori[index_pred]
  df=length(which(eta>0.5))
  BIC=log(error)+df*log(n)/n
  beta_interaction=mean_w_red[c(within_interaction_location_int,c((pp_int+1):b_s_int))]
  beta_main=mean_w_red[node_location_int]
  temp_2=list(mean_w_red=mean_w_red,beta_interaction=beta_interaction,beta_main=beta_main,BIC=BIC,df=df,eta_ori=eta_ori,gamma=gamma,index_path=index_path_ori,error=error)
  result_temp<-temp_2
  beta_interaction=mean_w_red[c(within_interaction_location_int,c((pp_int+1):b_s_int))]
  beta_main=mean_w_red[node_location_int]
  cat('df=',df,'\n')
  return(result_temp)

}

temp_function<-function(a,b,L){
  c=dim(L)[1]
  temp<-matrix(0,(c-a),(c-b))
  temp[(b-a),]=L[a,(b+1):c]
  temp1=diag(L[a,b],(c-b),(c-b))
  temp[(b-a+1):(c-a),]=temp1
  return(temp)
}

temp_function2<-function(L){
  s=1
  c=dim(L)[1]
  L_inter=matrix(0,c*(c-1)/2,c*(c-1)/2)
  for(j in 1:(c-1)){
    temp=NULL
    for(i in j:(c-1)){
      temp=cbind(temp,temp_function(j,i,L))
    }
    L_inter[s:(s+c-j-1),s:((c-1)*c/2)]=temp
    s=s+c-j
  }
  temp=L_inter+t(L_inter)
  for(i in 1:c*(c-1)/2){
    temp[i,i]=1
  }
  return(temp)
}






location<-function(p){
  if(p>1){
    temp=matrix(0,p,p)
    r=1
    for(i in 1:(p-1)){
      temp[i,(i+1):p]=c(r:(r+p-i-1))
      r=r+p-i
    }
  }else{
    temp=0
  }
  return(temp)
}






