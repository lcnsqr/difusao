typedef struct {
	int width, height, depth;
	int frameSize;
	char *frame;
	char *fgColor[2];
	char *bgColor;
} Video;

void videoBuild(Video* video, const int width, const int height, const int depth);
void videoGraphVector(Video* video, const float min, const float max, const float* v, const int m);

void videoBuild(Video* video, const int width, const int height, const int depth){
	video->width = width;
	video->height = height;
	video->depth = depth;
	video->frameSize = video->width*video->height*video->depth;

	video->frame = (char *)malloc(sizeof(char)*video->frameSize);
	video->fgColor[0] = (char *)malloc(sizeof(char)*video->depth);
	video->fgColor[1] = (char *)malloc(sizeof(char)*video->depth);
	video->bgColor = (char *)malloc(sizeof(char)*video->depth);

	video->bgColor[0] = 255;
	video->bgColor[1] = 255;
	video->bgColor[2] = 255;
	video->bgColor[3] = 255;
	video->fgColor[0][0] = 18;
	video->fgColor[0][1] = 10;
	video->fgColor[0][2] = 0;
	video->fgColor[0][3] = 255;
	video->fgColor[1][0] = 0;
	video->fgColor[1][1] = 100;
	video->fgColor[1][2] = 30;
	video->fgColor[1][3] = 255;
	#pragma omp parallel for
	for (int i = 0; i < video->frameSize; i += video->depth)
		for (int d = 0; d < video->depth; d++)
			video->frame[i+d] = video->bgColor[d];
}

void videoPoint(Video* video, const float min, const float max, const float* v){
	// Proporção entre domínios
	const float r = (float)video->width / fabs(max-min);;
	int x, y;
	x = floor(r*(v[0] - min));
	y = floor(r*(v[1] - min));
	// A orientação vertical é invertida na imagem
	y = video->height - y - 1;
	// Sair se ponto fora da imagem
	if ( x < 0 || x >= video->width || y < 0 || y >= video->height ) return;
	// Single pixel
	for (int d = 0; d < video->depth; d++)
		video->frame[video->depth * (y * video->width + x) + d] = video->fgColor[0][d];
	// Draw a square
	/*
	for (int sy = -1; sy < 2; sy++)
		for (int sx = -1; sx < 2; sx++)
			for (int d = 0; d < video->depth; d++)
				if ( y + sy > 0 && y + sy < video->height && x + sx > 0 && x + sx < video->width )
					video->frame[video->depth * ((y+sy) * video->width + (x+sx)) + d] = video->fgColor[0][d];
	*/
}

void videoPointVector(Video* video, const float min, const float max, const float* v, const int m){
	// Proporção entre domínios
	const float r = (float)video->width / fabs(max-min);;
	int x, y;
	#pragma omp parallel for private(x,y)
	for (int i = 0; i < 2*m; i += 2){
		x = floor(r*(v[i] - min));
		y = floor(r*(v[i+1] - min));
		// A orientação vertical é invertida na imagem
		y = video->height - y - 1;
		// Draw a square
		for (int sy = -1; sy < 1; sy++)
			for (int sx = -1; sx < 1; sx++)
				for (int d = 0; d < video->depth; d++)
					if ( y + sy > 0 && y + sy < video->height && x + sx > 0 && x + sx < video->width )
						video->frame[video->depth * ((y+sy) * video->width + (x+sx)) + d] = video->fgColor[0][d];
	}
}

void videoGraphVector(Video* video, const float min, const float max, const float* v, const int m){
	// Proporção entre a quantidade de pontos no vetor e a resolução da imagem
	double r = (double)m/video->width;
	int x, y;
	#pragma omp parallel for private(x,y)
	for (int i = 0; i < video->frameSize; i += video->depth){
		x = (i/video->depth) % video->width;
		y = (i/video->depth) / video->width;
		// A orientação vertical é invertida na imagem
		y = video->height - y - 1;
		// Cor de fundo
		for (int d = 0; d < video->depth; d++)
			video->frame[i+d] = video->bgColor[d];
		// Target
		/*
		if ( (int)floor((video->height-1)*(*video->target)) == y )
			for (int d = 0; d < video->depth; d++)
				video->frame[i+d] = video->fgColor[0][d];
		*/
		if ( (int) floor((video->height-1) * ( v[(int)floor(r*x)] - min ) / (max-min) ) == y )
			for (int d = 0; d < video->depth; d++)
				video->frame[i+d] = video->fgColor[1][d];
	}
}

void videoGraphVectorColor(Video* video, struct Cor* cor, const float min, const float max, const float* v, const int m){
	// Proporção entre a quantidade de pontos no vetor e a resolução da imagem
	const float r = (float)m/video->width;
	int x, y, q, f;
	char rgba[4];
	#pragma omp parallel for private(x,y,q,f,rgba)
	for (int i = 0; i < video->frameSize; i += video->depth){
		x = (i/video->depth) % video->width;
		y = (i/video->depth) / video->width;
		// A orientação vertical é invertida na imagem
		y = video->height - y - 1;
		// Valor relativo à posição x
		q = v[(int)floor(r*x)];
		// Valor da função no quadro
		f = (video->height-1) * ( (int)q - min ) / (max-min);
		// Cor de fundo
		/*
		for (int d = 0; d < video->depth; d++)
			video->frame[i+d] = video->bgColor[d];
		*/
		// Cor relativa ao valor da função no ponto atual
		corMap(cor, f, rgba);
		for (int d = 0; d < video->depth; d++)
			video->frame[i+d] = rgba[d];
		// Traço da função
		if ( f == y || f == y + 1 || f == y - 1 )
			for (int d = 0; d < video->depth; d++)
				video->frame[i+d] = 255 - rgba[d];
	}
}
