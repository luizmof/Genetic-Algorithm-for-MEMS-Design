%% Algoritmo de Sele��o Clonal - Otimiza��o de param�tros para modelagem de
%% MEMS aplicados a Energy Harvesting.
%% Trabalho Final - Professor Mesquita 2016.3
%% Aluno: Luiz Carlos Macedo de Oliveira Filho
%%DRE: 116217365

%% O programa ir� gerar valores para otimizar o projeto do MEMS para Energy Harvesting proposto inicialmente no meu projeto final e 
%%que ser� complementado em minha disserta��o de mestrado.
%%A parte das limita��es de valores para as vari�veisa serem
%%postas no ciclo de muta��o vem do conhecimento e heur�stica adquiridos pela Inicia��o
%%Cientif�ca e Projeto final com o Professor Moreir�o sobre o assunto, e a
%%proposta agora para o mestrado seria a fabrica��o de um modelo em
%%microeletr�nica, tendo em vista que um modelo de teste foi feito em impress�o 3D na escala de mm. Alguns parametros ser�o otimizados e futuramente
%%procura-se tamb�m explorar essa aplica��o para outros fatores do projeto.
%%A proposta ser� fazer com uma popula��o de 50 cromossomos de 49 bits selecionando
%%sempre os 10 melhores e avaliando eles segundo crit�rios estabelecidos por mim.

%% Defini��o do cromossomo a ser usado: [Ndentes | L | D | W | E | H | PHastes]
%% Ndentes -> N�mero de dentes das hastes, ser� modificado com a muta��o. Variando de 3 a 128, 2^7. Dos bits 1 a 7 da linha.
%% L -> Comprimento do dente, ser� modificado com a muta��o. Variando de 5 micro a 128 micro.Dos bits 8 a 14 da linha.
%% D -> Dist�ncia entre os dentes, ser� modificado com a muta��o. Variando de 6 a 128 micro. Dos bits 15 a 21 da linha.
%% W -> Largura dos dentes, ser� modificado com a muta��o. Variando de 6 a 128 micro. Dos bits 22 a 28 da linha.
%% E -> Largura das hastes, ser� modificado com a muta��o. Variando entre 6 a 128 micro. Dos bits 29 a 35 da linha.
%% H -> Altura da estrutura, ser� modificado com a muta��o. Variando entre 1 a 128 micro. Dos bits 36 a 42 da linha.
%% Phastes -> N�mero de hastes, ser� modificado com a muta��o. Variando de 3 a 128. Dos bits 43 a 49 da linha.
%% Cmax -> Capacit�ncia m�xima, n�o pode ser modificado com muta��o, pois � calculado. Cmax = (Eo*Smax*H*Phastes)/(D).
    %%Valor em torno de 150pF
%% Cmin -> Capacit�ncia m�nima, n�o pode ser modificado com muta��o, pois � calculado. Cmin = (Eo*Smin*H*Phastes)/(E) - Eo*(L+E)*H/E.
%% Smax -> Superf�cie de sobreposi��o para capacit�ncia m�xima, n�o pode ser modificado com muta��o, pois � calculado.
    %%Smax = (2*N-1)*(L-D)+2*(N-1)*W+(2*N-1)*alfa*D.
%% Smin -> Superf�cie de sobreposi��o para capacit�ncia minima, n�o pode ser modificado com muta��o, pois � calculado. 
    %%Smin = (N-1)*2*(D+W)-D+2*(L+E)+2*alfa*E.
    
%% Razao -> Razao entre as capacit�ncias Cmax e Cmin, n�o pode ser modificado com muta��o, pois � calculado. Razao = Cmax/Cmin. 
    %%Razao tem que ser maior que 2, objetivo ser� posto em 3 e quanto
    %%maior melhor, assim n�o tendo limite m�ximo. Ser� armazenado em um
    %%vetor separado de double.
    
 %% Taxa de muta��o em 60% dos genes dos cromossomos.
 %% Erro na raz�o das capacit�ncia m�ximas e m�nimas em 1%.
    
%--Primeiro Passo: Gerar a popula��o aleat�rio.--
erro = 100; %aApenas um valor inicial maior que o objetivo
alfa = 0.5;
Eo = 8.854*10^-12;
objetivo = 3.3; %Valor objetivo do erro desejado.
iteracao = 1;

binAux = zeros(1,5);
binAux2 = zeros(1,5);
binAux3 = zeros(1,4);
binAux4 = zeros(1,3);
binAux5 = zeros(1,2);
binAux6 = zeros(1,3);
binAux7 = zeros(1,4);


cromossomo = randi([0 1],50,49);
Razao = zeros(50,1);

while erro > 0.01 %Verifica��o do erro no loop principal.
%Pegando os valores da vari�veis 
    for indice = 1:50                
        %Os ifs a seguir servem para verificar se um dos valores � zero,
        %como tamb�m pode servir para setar valores minimos caso seja, como � o caso
        %da vari�vel D.
        Ndentes = bi2de(cromossomo(indice,[1:7]));
        if (Ndentes <= 1)            
            binAux = de2bi(3,7);
           	cromossomo(indice,[1:7]) = binAux;
            Ndentes = 3;
        end;
        
        L = bi2de(cromossomo(indice,[8:14]));
        if (L <= 0)            
            binAux2 = de2bi(5,7);
           	cromossomo(indice,[8:14]) = binAux2;
            L = 5;            
        end;
        
        D = bi2de(cromossomo(indice,[15:21]));        
        if (D < 6)           
           binAux3 = de2bi(6,7);
           cromossomo(indice,[15:21]) = binAux3;
           D = 6;
        end;
        
        W = bi2de(cromossomo(indice,[22:28]));
        if (W <= 0)
           binAux4 = de2bi(6,7);
           cromossomo(indice,[22:28]) = binAux4;
           W = 6;
        end;
        
        E = bi2de(cromossomo(indice,[29:35]));
        if (E <= 0)           
           binAux5 = de2bi(6,7);
           cromossomo(indice,[29:35]) = binAux5;
           E = 6;
        end;
        
        H = bi2de(cromossomo(indice,[36:42]));
        if (H <= 0)            
            binAux6 = de2bi(1,7);
            cromossomo(indice,[36:42]) = binAux6;
            H = 1;            
        end;
        
        Phastes = bi2de(cromossomo(indice,[43:49]));
        if (Phastes <= 1)            
            binAux7 = de2bi(3,7);
            cromossomo(indice,[43:49]) = binAux7;
            Phastes = 3;
        end;

        %Calculando os parametros fixos e de controle.

        SmaxAux = (2*Ndentes-1)*(L-D)*10^-6+2*(Ndentes-1)*W*10^-6+(2*Ndentes-1)*alfa*D*10^-6;
        SminAux = (Ndentes-1)*2*(D+W)*10^-6-D*10^-6+2*(L+E)*10^-6+2*alfa*E*10^-6;

        CmaxAux = (Eo*SmaxAux*H*(10^-6)*Phastes)/(D*(10^-6));
        CminAux = (Eo*SminAux*H*(10^-6)*Phastes)/(E*10^-6) - Eo*(L+E)*(10^-6)*H*(10^-6)/(E*10^-6);        
        
        RazaoAux = CmaxAux/CminAux;

        Razao(indice) = RazaoAux;        

    end;

%%-- Segundo Passo: Calculo do Erro e Sele��o dos elementos com as 10 melhores raz�es.
    
    RazaoSort = zeros(10,1);
    RazaoAux = sort(Razao,'descend');
    contadorAux = 1;
    for n = 1:50        
        if (isnan(RazaoAux(n)))
            contadorAux = contadorAux + 1;
        else
            RazaoSort(n) = RazaoAux(n);
        end;
    end;    
    selecao = 1;
    contador1 = 1;
    
    %Calculo do erro:
    
    dezMelhores = zeros(50,1);
    selecionados = zeros(10,49);
    while selecao<11
        %Verifica��o quanto a quantidade de elementos selecionados, pode ser que 10 n�o cheguem a atender os requisitos, pode ser que nenhum chegue. 
        if contador1 > 50
            naoZero = find(selecionados);
            localNaoZero = size(naoZero);
            tamanho = localNaoZero(1) + 1;
            if tamanho > 1
                while tamanho < 11   
                    selecionados(tamanho) = selecionados(1);
                    tamanho = tamanho + 1;                
                end;
            else
                % Caso nenhum tenha sido selecionado, ser� sorteado 10
                % outros para serem clonados e mutados, assim gerando uma
                % popula��o toda nova.
                for n=1:10
                    selecionados = randi([0 1],10,49);                    
                end;                                
            end;
            break;
        end;
        
        %Inicio da avalia��o e calculo de erro.        
        
        dezMelhores(contador1) = find(Razao == RazaoSort(contador1),1);                
        
        Ndentes = bi2de(cromossomo(dezMelhores(contador1),[1:7]));
        L = bi2de(cromossomo(dezMelhores(contador1),[8:14]));
        D = bi2de(cromossomo(dezMelhores(contador1),[15:21]));
        W = bi2de(cromossomo(dezMelhores(contador1),[22:28]));
        E = bi2de(cromossomo(dezMelhores(contador1),[29:35]));
        H = bi2de(cromossomo(dezMelhores(contador1),[36:42]));
        Phastes = bi2de(cromossomo(dezMelhores(contador1),[43:49]));
                                
        SmaxAux = (2*Ndentes-1)*(L-D)*10^-6+2*(Ndentes-1)*W*10^-6+(2*Ndentes-1)*alfa*D*10^-6;
        SminAux = (Ndentes-1)*2*(D+W)*10^-6-D*10^-6+2*(L+E)*10^-6+2*alfa*E*10^-6;
        
        testeCapacitorMax = (Eo*SmaxAux*H*(10^-6)*Phastes)/(D*(10^-6));         
        testeCapacitorMin = (Eo*SminAux*H*(10^-6)*Phastes)/(E*10^-6) - Eo*(L+E)*(10^-6)*H*(10^-6)/(E*10^-6);        
        %Inicio das verifica��es de sa�da do programa.
        if testeCapacitorMin > 0 %Devido a f�rmula  para Smin, pode ocorrer, o que seria uma dimens�o n�o v�lida.
            if testeCapacitorMax > 150*10^-12
                %Sele��o dos 10 melhores                
                selecionados(selecao) = cromossomo(dezMelhores(contador1));             
                selecao = selecao + 1;
                erro = abs((objetivo-RazaoSort(contador1))/(objetivo));
                if erro < 0.01
                    %Verifica��o do erro, pela minha proposta eu avalio
                    %se a raz�o est� num valor em torno de 3.3 com erro m�ximo de 1 % e a Cmax
                    %maior que 150pF.
                    theChosenOne = cromossomo(dezMelhores(contador1),1:49); % Solu��o encontrada.
                    break;
                end;
            end;
        end;
        contador1 = contador1 + 1;
    end;        
    
%%-- Terceiro Passo: Clonar os 10 melhores elementos
    
    contador2 = 1;
    for clonar = 1:50
        cromossomo(clonar) = selecionados(contador2);
        if contador2 == 10
            contador2 = 1;
        else
            contador2 = contador2 + 1;
        end;        
    end;
    
%%-- Quarto Passo: Muta��o nos genes dos cromossomos.
    sorteio = randi([1 49],30,1);
    
    for mutacao = 1:50
        %Selecionar os genes a serem alterados
        for indice2 = 1:30
            if (cromossomo(mutacao,sorteio(indice2)) == 1)
                cromossomo(mutacao,sorteio(indice2)) = 0;
            else
                cromossomo(mutacao,sorteio(indice2)) = 1;
            end;                    
        end;
    end;

iteracao = iteracao+1;
    
end;

%%-- Impress�o do Resultado:
Ndentes = bi2de(theChosenOne(1,[1:7]));
L = bi2de(theChosenOne(1,[8:14]))*10^-6;
D = bi2de(theChosenOne(1,[15:21]))*10^-6;
W = bi2de(theChosenOne(1,[22:28]))*10^-6;
E = bi2de(theChosenOne(1,[29:35]))*10^-6;
H = bi2de(theChosenOne(1,[36:42]))*10^-6;
Phastes = bi2de(theChosenOne(1,[43:49]));

%Calculando os parametros

SmaxAux = (2*Ndentes-1)*(L-D)+2*(Ndentes-1)*W+(2*Ndentes-1)*alfa*D;
SminAux = (Ndentes-1)*2*(D+W)-D+2*(L+E)+2*alfa*E;

CmaxAux = (Eo*SmaxAux*H*Phastes)/(D);
CminAux = (Eo*SminAux*H*Phastes)/(E) - Eo*(L+E)*H/(E);

RazaoAux = CmaxAux/CminAux;

Ndentes
L
D
W
E
H
Phastes
SmaxAux
SminAux
CmaxAux
CminAux
RazaoAux
erro