#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "SDL.h"
#include "SDL_ttf.h"
#include "matrix.h"
#include "color.h"
#include "video.h"

// Solução numérica da equação do calor pelo
// método Crank-Nicolson de diferenças finitas

// Extremidade inicial e final
#define POS0 0
#define POS1 100

// Menor e maior valores de leitura
#define MINVAL 0
#define MAXVAL 1

// Constantes de exibição para o modo interativo
#define DEPTH 4
#define WIDTH 800
#define HEIGHT 320

// Funções auxiliares para preenchimento de matriz
float zeros(int i, int j){
	return 0;
}
float eye(int i, int j){
	return ( i == j ) ? 1 : 0;
}

// Escrever vetor em arquivo
int saveArray(float* res, const char* file, const int m){
	int i = 0;
	FILE* f = fopen(file, "w");
	while ( i < m ){
		fprintf(f, "%.25f\n", res[i]);
		i++;
	}
	fclose(f);
	return i;
}

// Construção da matriz A em Aw_j+1 = Bw_j com 
// a condição de Dirichlet em a e Neuman em b.
// A: Matriz resultante
// lambda: constante de diferenças finitas
void leftA(Matrix* A, const float lambda){
	A->_[0] = 1;
	for (int i = 1; i < A->rows - 2; i++){
		A->_[i*A->cols+i-1] = -lambda/2;
		A->_[i*A->cols+i] = 1+lambda;
		A->_[i*A->cols+i+1] = -lambda/2;
	}
	A->_[(A->rows-2)*A->cols+A->cols-3] = -lambda/2;
	A->_[(A->rows-2)*A->cols+A->cols-2] = 1+lambda/2;
	A->_[(A->rows-2)*A->cols+A->cols-1] = -lambda/2;
	A->_[A->size-1] = 1;
}

// Construção da matriz B em Aw_j+1 = Bw_j com 
// a condição de Dirichlet em a e Neuman em b.
// B: Matriz resultante
// lambda: constante de diferenças finitas
void rightB(Matrix* B, const float lambda){
	B->_[0] = 1;
	for (int i = 1; i < B->rows - 2; i++){
		B->_[i*B->cols+i-1] = lambda/2;
		B->_[i*B->cols+i] = 1-lambda;
		B->_[i*B->cols+i+1] = lambda/2;
	}
	B->_[(B->rows-2)*B->cols+B->cols-3] = lambda/2;
	B->_[(B->rows-2)*B->cols+B->cols-2] = 1-lambda/2;
	B->_[(B->rows-2)*B->cols+B->cols-1] = lambda/2;
	B->_[B->size-1] = 1;
}

// Gerar inversa de matriz tridiagonal
// Res: Matriz resultante (inversa de A)
// A: Matriz tridiagonal
void invTri(Matrix* res, Matrix* A){
	// Matriz res deve ser identidade
	mtrxRebuildWith(res, A->rows, A->cols, &eye);
	// Copiar matriz, inversão altera seu conteúdo
	Matrix M;
	mtrxBuildNull(&M, A->rows, A->cols);
	mtrxEqual(&M, A);
	// Vetores auxiliares
	Matrix row[3], rowInv[3];
	mtrxBuild(&row[0], 1, M.cols);
	mtrxBuild(&row[1], 1, M.cols);
	mtrxBuild(&row[2], 1, M.cols);
	mtrxBuild(&rowInv[0], 1, M.cols);
	mtrxBuild(&rowInv[1], 1, M.cols);
	mtrxBuild(&rowInv[2], 1, M.cols);
	// Posição na matriz
	int p[2];
	p[0] = 0;
	p[1] = 0;
	// Matrix tridiagonal, com diagonal dominante estrita. 
	// Inversão por escalonamento simples.
	// Primeira linha já pronta
	for (int i = 1; i < M.rows - 1; i++){
		// i-ésima linha
		mtrxRow(&row[0], &M, i);
		// Na inversa
		mtrxRow(&rowInv[0], res, i);
		// Dividir pelo elemento em i-1
		mtrxScalar(&row[1], &row[0], 1.0/row[0]._[i-1]);
		// Na inversa
		mtrxScalar(&rowInv[1], &rowInv[0], 1.0/row[0]._[i-1]);
		// Subtrair com a linha anterior
		mtrxRow(&row[0], &M, i-1);
		mtrxMinus(&row[2], &row[1], &row[0]);
		// Na inversa
		mtrxRow(&rowInv[0], res, i-1);
		mtrxMinus(&rowInv[2], &rowInv[1], &rowInv[0]);
		// Dividir pelo elemento em i
		mtrxScalar(&row[0], &row[2], 1.0/row[2]._[i]);
		// Na inversa
		mtrxScalar(&rowInv[0], &rowInv[2], 1.0/row[2]._[i]);
		// Atualizar linha na Matriz
		p[0] = i;
		mtrxPaste(&M, &row[0], p);
		// Na inversa
		mtrxPaste(res, &rowInv[0], p);
	}
	// Finalizar inversão diagonal-superior resultante
	// Última linha já pronta
	for (int i = M.rows - 2; i > 0; i--){
		// i-ésima linha
		mtrxRow(&row[0], &M, i);
		// Na inversa
		mtrxRow(&rowInv[0], res, i);
		// (i+1)-ésima linha
		mtrxRow(&row[1], &M, i+1);
		// Na inversa
		mtrxRow(&rowInv[1], res, i+1);
		// Multiplicar pelo elemento acima na linha superior
		mtrxScalar(&row[2], &row[1], row[0]._[i+1]);
		// Na inversa
		mtrxScalar(&rowInv[2], &rowInv[1], row[0]._[i+1]);
		// Subtrair i-ésima com (i+1)-ésima
		mtrxMinus(&row[1], &row[0], &row[2]);
		// Na inversa
		mtrxMinus(&rowInv[1], &rowInv[0], &rowInv[2]);
		// Atualizar linha na Matriz
		p[0] = i;
		mtrxPaste(&M, &row[1], p);
		// Na inversa
		mtrxPaste(res, &rowInv[1], p);
	}
}

int main(int argc, char** argv){

	// Recolher opções e parâmetros
	int opt;
	// Constante da equação
	float alfa;
	// Estado depois de t segundos
	int t;
	if ( argc < 2 ){
		fprintf(stderr, "Utilização: %s [-c constante] [-t segundos] < inicial.txt > final.txt \n", argv[0]);
		exit(EXIT_FAILURE);
	}
	while ((opt = getopt(argc, argv, "c:t:")) != -1) {
		switch (opt) {
			case 'c':
				alfa = fabs(atof(optarg));
			break;
			case 't':
				t = abs(atoi(optarg));
			break;
			default: /* '?' */
				fprintf(stderr, "Utilização: %s [-c constante] [-t segundos] < inicial.txt > final.txt \n", argv[0]);
				exit(EXIT_FAILURE);
		}
	}

	// Extremidade inicial
	float a = POS0;
	// Extremidade final
	float b = POS1;
	// Valor minimo de leitura
	float min = MINVAL;
	// Valor máximo de leitura
	float max = MAXVAL;

	// Vetores antes e depois
	// Condição de Dirichlet no extremo a
	// Valor no extremo a definido pelo primeiro valor do vetor
	// Condição de Neuman no extremo b
	// Valor da derivada em relação a x no extremo b definido pelo último valor do vetor 
	Matrix w[2];
	mtrxReadFromSTDIN(&w[0]);
	// Total de pontos entre a e b (inclusivo)
	int m = w[0].rows;
	// Aumentar o vetor para incluir o valor 
	// da condição de Neuman no extremo b 
	mtrxRebuild(&w[0], m+1, 1);
	w[0]._[m] = 0;
	// Copiar valores de antes para depois
	mtrxBuild(&w[1], m+1, 1);
	mtrxEqual(&w[1], &w[0]);

	// Constante da equação
	if ( alfa == 0 ){
		fprintf(stderr, "A constante da equação não pode ser nula. \n");
		exit(EXIT_FAILURE);
	}
	// Espaçamento entre pontos 
	float h = (b-a)/m;
	// Tamanho do passo no tempo calculado a partir 
	// da constante alfa e do espaçamento h
	float k = pow(h,2)/alfa;
	// Limitar o passo no tempo em 1 segundo
	k = ( k > 1 ) ? 1 : k;
	// Constante auxiliar lambda
	float lambda = k*alfa/pow(h,2);
	//fprintf(stderr, "h = %.9f\nk = %.9f\n", h, k);

	// Matrizes da operação
	Matrix A, invA, B;
	// Iniciá-las com zeros para posterior preenchimento
	mtrxBuildWith(&A, m+1, m+1, &zeros);
	mtrxBuildWith(&invA, m+1, m+1, &eye);
	mtrxBuildWith(&B, m+1, m+1, &zeros);
	// Preenchimento das matrizes tridiagonais
	leftA(&A, lambda);
	rightB(&B, lambda);
	// Inverter A
	invTri(&invA, &A);
	// Construir a matriz inv(A)*B
	Matrix C;
	mtrxBuild(&C, m+1, m+1);
	mtrxMul(&C, &invA, &B);
	// Matrizes A, B, inv(A) não mais necessárias
	mtrxDiscard(&A);
	mtrxDiscard(&invA);
	mtrxDiscard(&B);

	// Tempo acumulado das iterações
	float ts = 0;
	// Se t > 0, modo não interativo
	if ( t > 0 ){
		// Contar os passos no tempo até atingir o instante desejado
		while (ts <= t){
			// Próximo resultado a partir do estado atual
			mtrxMul(&w[1], &C, &w[0]);
			mtrxEqual(&w[0], &w[1]);
			ts += k;
		}
		// Vetor resultante na saída padrão
		for (int i = 0; i < m; i++) fprintf(stdout, "%.25f\n", w[0]._[i]);
		// Encerrar
		return 0;
	}

	// Modo gráfico interativo
	int width = WIDTH;
	int height = HEIGHT;
	// Estrutura que gera e armazena os quadros da animação
	Video video;
	videoBuild(&video, width, height, DEPTH);

	// Exibição e interação por SDL
	SDL_Window *window;
	SDL_Renderer *renderer;
	SDL_Texture *texture;
	SDL_Event event;

	if (SDL_Init(SDL_INIT_VIDEO) < 0) {
		SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Erro ao iniciar o SDL: %s", SDL_GetError());
		exit(-1);
	}

	if (SDL_CreateWindowAndRenderer(width, height, SDL_WINDOW_RESIZABLE, &window, &renderer)) {
		SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Impossível criar a janela e o renderizador: %s", SDL_GetError());
		exit(-1);
	}
	SDL_SetWindowTitle(window, "Difusão por Diferenças Finitas");

	texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STREAMING, width, height);

	// Renderizar texto (tempo decorrido)
	if ( TTF_Init() < 0 ){
		SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Erro ao executar TTF_Init: %s", TTF_GetError());
	}
	TTF_Font* txtFont = TTF_OpenFont("FreeSans.ttf", 18);
	SDL_Color txtColorF = {255, 255, 255, 255};
	SDL_Color txtColorB = {54, 78, 89, 255};
	char *txt = malloc(256);;
	sprintf(txt, "%.1f", 0);
	SDL_Surface* txtSurface = TTF_RenderText_Shaded(txtFont, txt, txtColorF, txtColorB);
	SDL_Texture* txtTexture = SDL_CreateTextureFromSurface(renderer, txtSurface);
	int txtWidth, txtHeight;
	SDL_QueryTexture(txtTexture, NULL, NULL, &txtWidth, &txtHeight);
	SDL_Rect txtRect;
	txtRect.x = width - txtWidth;
	txtRect.y = height - txtHeight;
	txtRect.w = txtWidth;
	txtRect.h = txtHeight;

	// Estrutura de mapeamento posição -> cor
	struct Cor cor;
	// Total de cores equivale à altura da imagem
	corBuild(&cor, height);

	// Posição do cursor relativo à janela
	int curPos[2];
	// Loop de exibição
	while (1){
		// Gerar gráfico
		videoGraphVectorColor(&video, &cor, min, max, w[0]._, m);

		// Atualizar exibição
		SDL_UpdateTexture(texture, NULL, video.frame, DEPTH * width * sizeof(char));
		SDL_RenderClear(renderer);
		SDL_RenderCopy(renderer, texture, NULL, NULL);

		// Atualizar exibição do tempo
		sprintf(txt, "Tempo decorrido: %.1f segundos", ts);
		SDL_FreeSurface(txtSurface);
		txtSurface = TTF_RenderText_Shaded(txtFont, txt, txtColorF, txtColorB);
		SDL_DestroyTexture(txtTexture);
		txtTexture = SDL_CreateTextureFromSurface(renderer, txtSurface);
		SDL_QueryTexture(txtTexture, NULL, NULL, &txtWidth, &txtHeight);
		txtRect.x = width - txtWidth;
		txtRect.y = height - txtHeight;
		txtRect.w = txtWidth;
		txtRect.h = txtHeight;
		SDL_RenderCopy(renderer, txtTexture, NULL, &txtRect);

		SDL_RenderPresent(renderer);

		// Recolher eventos
		SDL_PollEvent(&event);
		if (event.type == SDL_QUIT){
			break;
		}
		else if( event.type == SDL_KEYDOWN ){
			if ( event.key.keysym.sym == SDLK_q ){ 
				// Tecla "q" encerra
				break;
			}
			else if ( event.key.keysym.sym == SDLK_UP && w[0]._[0] < max ){ 
				// Seta pra cima aumenta o valor no extremo a
				w[0]._[0] += (max-min)/height;
			}
			else if ( event.key.keysym.sym == SDLK_DOWN && w[0]._[0] > min ){ 
				// Seta pra baixo reduz o valor no extremo a
				w[0]._[0] -= (max-min)/height;
			}
		}
		SDL_PumpEvents();
		if (SDL_GetMouseState(&curPos[0], &curPos[1]) & SDL_BUTTON(SDL_BUTTON_LEFT)) {
			// Se botão do mouse pressionado, aproximar o valor 
			// no extremo para a posição vertical do cursor
			w[0]._[0] += ((1 - (float)curPos[1]/height) * (max-min) - w[0]._[0]);
		}

		// Próximo resultado a partir do estado atual
		mtrxMul(&w[1], &C, &w[0]);
		mtrxEqual(&w[0], &w[1]);
		w[0]._[m] = 0;

		// Tempo decorrido
		ts += k;
	}

	// Encerrar SDL
	TTF_CloseFont(txtFont);
	SDL_DestroyTexture(texture);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);

	SDL_Quit();

	// Salvar estado
	//saveArray(w[0]._, argv[1], m);

	// Exibir vetor resultante na saída padrão
	for (int i = 0; i < m; i++) fprintf(stdout, "%.25f\n", w[0]._[i]);

	return 0;
}
