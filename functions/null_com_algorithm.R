null_com_algorithm <- function(com){
  
  #beta_dev_testset.seed(seed)
  
  #gerando matriz de presença e ausencia
  com_pa <- decostand(com, method = "pa")
  
  #Gerando vetor com numero de ocorrencias de cada especie
  ocorrences <- colSums(com_pa)
  
  #Gerando nova matriz
  null_com_pa <- matrix(data = 0, nrow = nrow(com_pa), ncol = ncol(com_pa))
  
  #Sorteando ocorrencias em cada coluna 
  lines <- 1:nrow(null_com_pa)#Criando objeto que guarda identidade das linhas (comunidades)
  for(i in 1:ncol(null_com_pa)){ #For loop para cada coluna
    if(ocorrences[i] > 0){ #Se número de ocorrências da espécie for maior que zero...
      new_presences <- sample(lines, ocorrences[i], replace = FALSE) #Sorteia novas linhas para cada presença
      null_com_pa[new_presences,i] <- 1 #Adiciona presença nas linhas sorteadas
    }
  }
  
  
  #Matrix com presenças embaralhadas pronta.
  
  #Gerando vetor com tamanhos das comunidades
  com_size <- rowSums(com)
  
  #Gerando vetor com tamanhos das populações
  pop_size <- colSums(com)
  
  #N total de individuos
  total_size <- sum(pop_size)
  
  #probabilidades para cada especie
  p_sp <- pop_size/total_size
  
  #Geranto vetor com nova riqueza da comunidae nula
  richness <- rowSums(null_com_pa)
  
  #Reduzindo o numero de especies do número de individuos para garantir ao menos um individuo em cada presença da matriz
  new_com_size <- com_size - richness
  
  
  ###########################333
  sum(com_pa) == sum(null_com_pa)# Daqui pra tras tudo ok
  ##############################
  
  
  #Agora é possível que algumas linhas tenham menos individuos do que presenças, fazendo com que o new_com_size seja negativo
  if(any(new_com_size < 0)){ #Algum tamanho de comunidade negativo? Se sim, segue
    
    negatives <- which(new_com_size < 0) #Identifica quais comunidades tem tamanho (subtraído da riqueza) menor que zero
    n_negatives <- length(negatives) #Estabelece o número de comunidades com tamanho negativo
    
    while(n_negatives > 0){ #while loop que se repete até não ter nenhuma comunidade com tamanho negativo
      
      #negatives <- which(new_com_size < 0)
      #n_negatives <- length(negatives)
      
      
     
      
      #Agora precisamos, em cada uma das linhas, retirar presenças e colocar em outras linhas até que tenhamos pelo menos um numero de presença igual ao de individuos
      
      for(i in 1:length(negatives)){ #for loop que se repete para o número de LINHAS com tamanho negativo. i = LINHA
        presences <- null_com_pa[negatives[i],] #Cria vetor com as presenças ou ausencias da comunidade com valores de tamanho negativo
        sps <- which(presences!= 0) #cria vetor com as localizações nas linhas das presenças 
        sps <- sps[colSums(null_com_pa[,sps]) < nrow(null_com_pa)] #Remover especies com presenças saturadas (espécies presentes em todas as comunidades, não tem como adicionar presença para essas)
        abs_neg <- abs(new_com_size[negatives[i]]) #Cria objeto com o NÚMERO DE PRESENÇAS que deverão ser adicionadas na comunidade com tamanho negativo, e numero de ausências adicionadas a outras comunidades
        if(length(sps) < 2){
          remove <- sps
        }else{
          remove <- sample(sps, abs_neg) #Sorteia localização (COLUNAS onde a espécies será retirada)
        }
        
        for(j in 1:abs_neg){# for loop para cada presença que será adicionada ou retirada
          sp_absences <- which(null_com_pa[,remove[j]] == 0) #objeto com as identidades das LINHAS onde espécie está ausente
          if(length(sp_absences) < 2){
            sampled_sp_absences <- sp_absences
          }else{
            sampled_sp_absences <- sample(sp_absences, 1) 
          }          
          null_com_pa[,remove[j]][sampled_sp_absences] <- 1 #Adicionando presença na linha sorteada, 
        }
        null_com_pa[negatives[i], remove] <- 0 #Remove onde tinha
      }
      
      richness <- rowSums(null_com_pa) #Nova riqueza após remover os negativos
      new_com_size <- com_size - richness #Novo tamahho de comunidade após remover as riquezas
      negatives <- which(new_com_size < 0) #Quais são os negativos após remover anteriormente
      n_negatives <- length(negatives) #Número de negativos após remover anteriormente
    }
    
    
  }
  
  
  
  
  #Em alguns casos é possível que o modelo nulo deixe linhas sem presenças, precisamos corrigir isso 
  while(any(richness  == 0)){
    problem <- which(richness == 0) #identifica quais linhas tem riqueza igual a zero
    
    #Agora precisamos, em cada uma das linhas, retirar presenças e colocar em outras linhas até que tenhamos pelo menos um numero de presença igual ao de individuos
    for(i in 1:length(problem)){ #loop para numero de linhas com riqueza igual a zero
      
      
      ###### PROBLEMA AQUI #################
      #####################################
      
      #Adicionando nova presença em espécies que já existem
      ocorrences_no_zero <- which(ocorrences != 0) #identificando espécies que ja existem
      if(length(ocorrences_no_zero) < 2){
        new_presence <- ocorrences_no_zero
      }else{
        new_presence <- sample(ocorrences_no_zero, 1) #sorteando uma para receber a nova presença
      }
      
      #retirando uma da espécie onde ela foi adicionada
      occurrences_sp_interest <- which(null_com_pa[,new_presence] == 1) #gerando vetor com os locais das presenças da espécie a qual as ocorrencias serão infladas
      
      if(length(occurrences_sp_interest) < 2){
        remove <- occurrences_sp_interest
      }else{
        remove <- sample(c(occurrences_sp_interest), size = 1) #sorteando um local para ela deixar de existir
      }
      
      
      null_com_pa[remove,new_presence] <- 0 #removendo a presença
      null_com_pa[problem[i], new_presence] <- 1 #recebendo a nova presença
      
      null_com_pa[,new_presence]
      #####################################
    }
    
    #Recapturando a riqueza
    richness <- rowSums(null_com_pa)
    
    #Reduzindo o numero de especies do número de individuos para garantir ao menos um individuo em cada presença da matriz
    new_com_size <- com_size - richness
    
    
    

    
    
    #Precisamos rodar de novo todo o código que verifica a presença de tamanhos negativos de comunidades
    if(any(new_com_size < 0)){ #Algum tamanho de comunidade negativo? Se sim, segue
      
      negatives <- which(new_com_size < 0) #Identifica quais comunidades tem tamanho (subtraído da riqueza) menor que zero
      n_negatives <- length(negatives) #Estabelece o número de comunidades com tamanho negativo
      
      while(n_negatives > 0){ #while loop que se repete até não ter nenhuma comunidade com tamanho negativo
        
        #negatives <- which(new_com_size < 0)
        #n_negatives <- length(negatives)
        
        
        
        
        #Agora precisamos, em cada uma das linhas, retirar presenças e colocar em outras linhas até que tenhamos pelo menos um numero de presença igual ao de individuos
        
        for(i in 1:length(negatives)){ #for loop que se repete para o número de LINHAS com tamanho negativo. i = LINHA
          presences <- null_com_pa[negatives[i],] #Cria vetor com as presenças ou ausencias da comunidade com valores de tamanho negativo
          sps <- which(presences!= 0) #cria vetor com as localizações nas linhas das presenças 
          sps <- sps[colSums(null_com_pa[,sps]) < nrow(null_com_pa)] #Remover especies com presenças saturadas (espécies presentes em todas as comunidades, não tem como adicionar presença para essas)
          abs_neg <- abs(new_com_size[negatives[i]]) #Cria objeto com o NÚMERO DE PRESENÇAS que deverão ser adicionadas na comunidade com tamanho negativo, e numero de ausências adicionadas a outras comunidades
          if(length(sps) < 2){
            remove <- sps
          }else{
            remove <- sample(sps, abs_neg) #Sorteia localização (COLUNAS onde a espécies será retirada)
          }
          
          for(j in 1:abs_neg){# for loop para cada presença que será adicionada ou retirada
            sp_absences <- which(null_com_pa[,remove[j]] == 0) #objeto com as identidades das LINHAS onde espécie está ausente
            if(length(sp_absences) < 2){
              sampled_sp_absences <- sp_absences
            }else{
              sampled_sp_absences <- sample(sp_absences, 1) 
            }          
            null_com_pa[,remove[j]][sampled_sp_absences] <- 1 #Adicionando presença na linha sorteada, 
          }
          null_com_pa[negatives[i], remove] <- 0 #Remove onde tinha
        }
        
        richness <- rowSums(null_com_pa) #Nova riqueza após remover os negativos
        new_com_size <- com_size - richness #Novo tamahho de comunidade após remover as riquezas
        negatives <- which(new_com_size < 0) #Quais são os negativos após remover anteriormente
        n_negatives <- length(negatives) #Número de negativos após remover anteriormente
      }
      
      
    }
    
  }
  

  
  #Distribuindo os individuos nas presenças
  null_com <- null_com_pa
  for(i in 1:nrow(null_com)){
    cols_with_presence <- which(null_com[i,] != 0)
    cols <- 1:ncol(null_com)
    if(length(cols_with_presence) == 1){
      sp_new_ind <- rep(cols_with_presence, new_com_size[i])
    }else{
      sp_new_ind <- sample(cols_with_presence, size = new_com_size[i], prob = p_sp[cols_with_presence], replace = TRUE)
    }
    
    for(sp in cols_with_presence){
      new_ab <- length(which(sp_new_ind == sp))
      null_com[i,sp] <- null_com[i,sp] + new_ab
    }
  }
  
  
  ###########################333
  sum(com_pa) == sum(null_com_pa)
  ##############################
  
  
  
  null_com <- data.frame(null_com)
  rownames(null_com) <- rownames(com)
  colnames(null_com) <- colnames(com)
  
  return(null_com)
  
}

