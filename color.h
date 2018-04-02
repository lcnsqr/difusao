// Mapeamento de número inteiro para cor
struct Cor { 
  // Numero de cores disponveis
  int estagios;
  // Pontos de passagem
  int verticesQuant;
  // Intervalos entre vertices
  int intervQuant;
  // Tamanho de cada intervalo
  float intervTaman;
  // Os vertices
  int** vertices;
  // Direcoes (vertices - 1)
  int** direcoes;
};

void corBuild(struct Cor* cor, int estagios){
  cor->estagios = estagios;
  cor->verticesQuant = 5;
  cor->intervQuant = cor->verticesQuant - 1;
  cor->intervTaman = cor->estagios / cor->intervQuant;

  cor->vertices = (int**) malloc(sizeof(int*) * cor->verticesQuant);
  for ( int i = 0; i < cor->verticesQuant; i++ ){
    cor->vertices[i] = (int*) malloc(sizeof(int) * 4);
  }

  cor->direcoes = (int**) malloc(sizeof(int*) * cor->intervQuant);
  for ( int i = 0; i < cor->intervQuant; i++ ){
    cor->direcoes[i] = (int*) malloc(sizeof(int) * 4);
  }

  // Pontos de passagem de cor
  cor->vertices[0][0] = 0;
  cor->vertices[0][1] = 0;
  cor->vertices[0][2] = 255;
  cor->vertices[0][3] = 255;

  cor->vertices[1][0] = 0;
  cor->vertices[1][1] = 255;
  cor->vertices[1][2] = 255;
  cor->vertices[1][3] = 255;

  cor->vertices[2][0] = 0;
  cor->vertices[2][1] = 255;
  cor->vertices[2][2] = 0;
  cor->vertices[2][3] = 255;

  cor->vertices[3][0] = 255;
  cor->vertices[3][1] = 255;
  cor->vertices[3][2] = 0;
  cor->vertices[3][3] = 255;

  cor->vertices[4][0] = 255;
  cor->vertices[4][1] = 0;
  cor->vertices[4][2] = 0;
  cor->vertices[4][3] = 255;

  // Definir direcoes entre vertices
  for ( int i = 0; i < cor->intervQuant; i++ ){
    for ( int k = 0; k < 4; k++ ){
      cor->direcoes[i][k] = cor->vertices[i + 1][k] - cor->vertices[i][k];
    }
  }
}

void corMap(struct Cor* cor, int c, char* rgba){
  int k, p;
  k = cor->intervQuant * c / cor->estagios;
  p = c - cor->intervTaman * k;
  for ( int j = 0; j < 4; j++ ){
    rgba[j] = cor->vertices[k][j] + p * cor->direcoes[k][j] / cor->intervTaman;
  }
}
