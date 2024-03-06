/* STIUJ Emanuel-Stefan - 311CC */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct pixel {
    unsigned char red;
    unsigned char green;
    unsigned char blue;
} pixel;

typedef struct imgNode {
    // structura unui nod apartinand arborelui
    // pixel == NULL pentru nodurile interne
    // topLeft == topRight == bottomRight == bottomLeft == NULL pentru frunze
    pixel *rgbValues;
    struct imgNode *topLeft;
    struct imgNode *topRight;
    struct imgNode *bottomRight;
    struct imgNode *bottomLeft;
    unsigned int depth;
} imgNode;

/* S-au implementat structurile si functiile pentru coada 
   cu ajutorul listelor simplu inlantuite. Acestea se utilizeaza la
   cerinta 2 la parcuregrea pe nivel a arborelui si la cerinta 3 la
   construirea arborelui in ordinea corecta */

typedef struct queueNode {
    struct queueNode *next;
    imgNode *nodeAddress;
} queueNode;

typedef struct Queue {
    queueNode *front;
    queueNode *rear;
    int size;
} Queue;

Queue* initQueue(void)
{
    Queue *newQueue = malloc(sizeof(Queue));
    newQueue->front = NULL;
    newQueue->rear = NULL;
    newQueue->size = 0;
    return newQueue;
}

int isQueueEmpty(Queue *q)
{	
	if (!(q->front)) {
		return 1;
	}
	return 0;
}

void enqueue(Queue *q, imgNode *nodeAddress)
{	
	queueNode *newNode = malloc(sizeof(queueNode));
	newNode->nodeAddress = nodeAddress;
	newNode->next = NULL;

	if (q->front == NULL) {
		q->rear = newNode;
		q->front = newNode;
		q->size++;
	}
	else {
		q->rear->next = newNode;
		q->rear = newNode;
		q->size++;
	}
}

imgNode* dequeue(Queue* q)
{	
	queueNode *p = q->front;
	q->front = q->front->next;
	q->size--;

	imgNode *nodeAddress = p->nodeAddress;
	free(p);

	if (!q->front) {
		q->rear = NULL;
		q->front = NULL;
	}
	return nodeAddress;
}

void destroyQueue(Queue *q)
{	
	queueNode *p = q->front;
	queueNode *t = NULL;

	while (p) {
		t = p;
		p = p->next;
		free(t);
	}
	free(q);	
}

pixel** mallocPixelMatrix(unsigned int size)
{
    // functia aloca memorie pentru matricea de pixeli si returneaza adresa ei
    pixel** newMatrix = malloc(size * sizeof(pixel *));
    int i;
    for (i = 0; i < size; i++) {
        newMatrix[i] = malloc(size * sizeof(pixel));
    }
    return newMatrix;
}

void freePixelMatrix(pixel **pixelMatrix, unsigned int size)
{
    // functia elibereaza memoria pentru o matrice de pixeli
    int i;
    for (i = 0; i < size; i++) {
        free(pixelMatrix[i]);
    }
    free(pixelMatrix);
}

void writePixel(pixel *myPixel, FILE *outputFile)
{
    // functia scrie in fisier datele unui pixel (red, green, blue)
    fwrite(&myPixel->red, 1, 1, outputFile);
    fwrite(&myPixel->green, 1, 1, outputFile);
    fwrite(&myPixel->blue, 1, 1, outputFile);
}

pixel readPixel(FILE *inputFile)
{
    // functia citeste intr-o structura 'pixel' datele dintr-un fisier
    pixel newPixel;
    fread(&newPixel.red, 1, 1, inputFile);
    fread(&newPixel.green, 1, 1, inputFile);
    fread(&newPixel.blue, 1, 1, inputFile);
    return newPixel;
}

void writeNodeData(imgNode *node, FILE *outputFile)
{
    // functia este utilizata la cerinta 2
    // aceasta scrie in fisier datele unui nod din arbore in functie de tipul acestuia
    unsigned char one = 1, zero = 0;
    if (node->rgbValues) { 
        fwrite(&one, 1, 1, outputFile);
        writePixel(node->rgbValues, outputFile);           
    } else {
        fwrite(&zero, 1, 1, outputFile);
    }
}

void freeTree(imgNode *node)
{
    // functia elibeareaza memoria pentru un arbore dat ca parametru
    if (!node) {
        return;
    }

    freeTree(node->topLeft);
    freeTree(node->topRight);
    freeTree(node->bottomRight);
    freeTree(node->bottomLeft);

    if (node->rgbValues) {
        free(node->rgbValues);
    }
    free(node);
}

void levelOrder(imgNode *node, Queue *myQueue, FILE *outputFile, unsigned int finalSize)
{
    // functia este folosita la cerinta 2
    // se parcurge arborele pe nivel cu ajutorul unei cozi 
    /* coada va retine adresele nodurilor ale caror date
    urmeaza sa fie scrise in fisier */

    unsigned char one = 1;
    unsigned char zero = 0;
    fwrite(&finalSize, sizeof(unsigned int), 1, outputFile);

    if (node->rgbValues) {
        fwrite(&one, 1, 1, outputFile);
        writePixel(node->rgbValues, outputFile);
    } else {
        fwrite(&zero, 1, 1, outputFile);
        enqueue(myQueue, node);
        while (!isQueueEmpty(myQueue)) {
            node = dequeue(myQueue);
            if (node->topLeft) {
                writeNodeData(node->topLeft, outputFile);
                enqueue(myQueue, node->topLeft);
            }
            if (node->topRight) {
                writeNodeData(node->topRight, outputFile);
                enqueue(myQueue, node->topRight);
            }
            if (node->bottomRight) {
                writeNodeData(node->bottomRight, outputFile);
                enqueue(myQueue, node->bottomRight);
            }
            if (node->bottomLeft) {
                writeNodeData(node->bottomLeft, outputFile);
                enqueue(myQueue, node->bottomLeft);
            }
        }
    }
}

pixel** readDataFromPPM(FILE *inputFile, unsigned int *size)
{
    /* functia citeste datele dintr-un fisier ppm si returneaza
    matricea de pixeli (imaginea) */

    int width, height, maxValue;
    char fileFormat[3];

    fscanf(inputFile, "%s", fileFormat);
    fscanf(inputFile, "%d%d", &width, &height);
    fscanf(inputFile, "%d", &maxValue);
    fseek(inputFile, 1, SEEK_CUR);

    pixel **myImage = mallocPixelMatrix(width);

    int i, j;
    for (i = 0; i < width; i++) {
        for (j = 0; j < width; j++) {
            myImage[i][j] = readPixel(inputFile);
        }
    }    
    *size = width;

    return myImage;
}

int pow2(int exp)
{
    // functia calculeaza 2 la puterea unui exponent primit ca parametru
    int res = 1;
    int i;
    for (i = 0; i < exp; i++) {
        res = res * 2;
    }
    return res;
}

int getLog2(int n)
{
    // functia calculeaza logaritm in baza 2 dintr-un numar dat ca parametru
    int p = 1;
    while (pow2(p) != n) {
        p++;
    }
    return p;
}

long long getMean(pixel **myImage, unsigned int x, unsigned int y, unsigned int size, pixel *averagePixel)
{
    // functia calculeaza valoarea 'mean' pe o submatrice si urmeaza sa fie returnata
    // variabila 'averagePixel' retine media culorilor red, green, blue din submatrice
    unsigned long long redAverage = 0, greenAverage = 0, blueAverage = 0;
    unsigned int i, j;
    for (i = x; i < x + size; i++) {
        for (j = y; j < y + size; j++) {
            redAverage = redAverage + myImage[i][j].red;
            greenAverage = greenAverage + myImage[i][j].green;
            blueAverage = blueAverage + myImage[i][j].blue;
        }
    }
    redAverage = redAverage / (size * size);
    greenAverage = greenAverage / (size * size);
    blueAverage = blueAverage / (size * size);
    
    long long mean = 0, redMean, greenMean, blueMean;
    for (i = x; i < x + size; i++) {
        for (j = y; j < y + size; j++) {
            redMean = (redAverage - myImage[i][j].red) * (redAverage - myImage[i][j].red);
            greenMean = (greenAverage - myImage[i][j].green) * (greenAverage - myImage[i][j].green);
            blueMean = (blueAverage - myImage[i][j].blue) * (blueAverage - myImage[i][j].blue);
            mean = mean + redMean + greenMean + blueMean;
        }
    }
    mean = mean / (3 * size * size);

    averagePixel->red = redAverage;
    averagePixel->green = greenAverage;
    averagePixel->blue = blueAverage;

    return mean;
}

imgNode* recursiveBuild(pixel **myImage, unsigned int size, unsigned int x, unsigned int y, unsigned int factor, unsigned int depth, 
                        unsigned int *maxDepth, unsigned int *numOfLeaves, unsigned int *maxBlockSize)
{
    // functia parcurge imaginea (matricea de pixeli) si construieste recursiv arborele
    if (depth < *maxDepth) {
        *maxDepth = depth;
    }

    imgNode *newNode = NULL;

    if (depth == 0) {
        // cazul in care se size == 1 adica submatricea este un pixel
        // se creaza un nod in arbore care va fi frunza si se va returna
        *numOfLeaves = *numOfLeaves + 1;

        newNode = malloc(sizeof(imgNode));
        newNode->rgbValues = malloc(sizeof(pixel));

        *(newNode->rgbValues) = myImage[x][y];
        newNode->topLeft = newNode->topRight = newNode->bottomRight = newNode->bottomLeft = NULL;

        return newNode;
    }

    pixel *averagePixel = malloc(sizeof(pixel));
    if (getMean(myImage, x, y, size, averagePixel) <= factor) {
        // cazul in care mean <= factor
        // se creaza un nou nod in arbore care va fi frunza
        // acesta va retine media culorilor din submatricea parcursa anterior
        if (*maxBlockSize < size) {
            *maxBlockSize = size;
        }
        *numOfLeaves = *numOfLeaves + 1;

        newNode = malloc(sizeof(imgNode));
        newNode->rgbValues = averagePixel;
        newNode->topLeft = newNode->topRight = newNode->bottomRight = newNode->bottomLeft = NULL;

        return newNode;
    } else {
        free(averagePixel);
    }

    newNode = malloc(sizeof(imgNode));
    newNode->rgbValues = NULL;

    /* in cazul in care nu se indeplineste niciuna dintre conditiile de mai sus se va
    diviza submatricea actuala in alte 4 blocuri, reluandu-se verificarile recursiv */
    newNode->topLeft = recursiveBuild(myImage, size / 2, x, y, factor, depth - 1, maxDepth, numOfLeaves, maxBlockSize);
    newNode->topRight = recursiveBuild(myImage, size / 2, x, y + size / 2, factor, depth - 1, maxDepth, numOfLeaves, maxBlockSize);
    newNode->bottomRight = recursiveBuild(myImage, size / 2, x + size / 2, y + size / 2, factor, depth - 1, maxDepth, numOfLeaves, maxBlockSize);
    newNode->bottomLeft = recursiveBuild(myImage, size / 2, x + size / 2, y, factor, depth - 1, maxDepth, numOfLeaves, maxBlockSize);
    return newNode;
}

void initialize(pixel **myImage, unsigned int size, imgNode *node, unsigned int x, unsigned int y)
{
    /* functia se va folosi la cerinta 3 la construirea imaginii decomprimate, 
    parcurgandu-se arborele rezultat in urma citirii fisierului comprimat */
    if (node == NULL) {
        return;
    }

    if (node->rgbValues) {
        int i, j;
        for (i = x; i < x + size; i++) {
            for (j = y; j < y + size; j++) {
                myImage[i][j] = *(node->rgbValues);
            }
        }
    }

    initialize(myImage, size / 2, node->topLeft, x, y);
    initialize(myImage, size / 2, node->topRight, x, y + size / 2);
    initialize(myImage, size / 2, node->bottomRight, x + size / 2, y + size / 2);
    initialize(myImage, size / 2, node->bottomLeft, x + size / 2, y);
}

void decomprimare(FILE *inputFile, FILE *outputFile)
{
    /* avand in vedere ca fisierul comprimat retine datele despre arbore in ordinea
    parcurgerii pe nivel, se va uitiliza o coada in construirea acestuia */

    Queue *myQueue = initQueue();
    
    imgNode *myTree = malloc(sizeof(imgNode));

    myTree->rgbValues = NULL;
    myTree->topLeft = myTree->topRight = myTree->bottomRight = myTree->bottomLeft = NULL;
    myTree->depth = 0;
    unsigned char nodeType;
    unsigned int size;

    fread(&size, sizeof(unsigned int), 1, inputFile);
    fread(&nodeType, 1, 1, inputFile);

    if (nodeType) {
        // cazul in care radacina este o frunza
        myTree->rgbValues = malloc(sizeof(pixel));
        *(myTree->rgbValues) = readPixel(inputFile);
        myTree->topLeft = myTree->topRight = myTree->bottomLeft = myTree->bottomRight = NULL;
        myTree->depth = 1;
    } else {
        // cazul in care radacina este nod intern
        enqueue(myQueue, myTree);

        while(!isQueueEmpty(myQueue)) {
            /* avand in vedere ca in coada se vor retine doar noduri interne se va aloca de fiecare data
            memorie pentru cei 4 descendenti */

            imgNode *tempNode = dequeue(myQueue);

            tempNode->topLeft = malloc(sizeof(imgNode));
            tempNode->topRight = malloc(sizeof(imgNode));
            tempNode->bottomRight = malloc(sizeof(imgNode));
            tempNode->bottomLeft = malloc(sizeof(imgNode));

            tempNode->topLeft->depth = tempNode->topRight->depth = tempNode->bottomRight->depth = tempNode->bottomLeft->depth = tempNode->depth + 1;

            /* in functie de tipul nodului care este citit din fisier, se vor efectua operatiile necesare in construirea corecta a arborelui*/

            fread(&nodeType, 1, 1, inputFile);
            if (nodeType) {
                tempNode->topLeft->rgbValues = malloc(sizeof(pixel));
                *(tempNode->topLeft->rgbValues) = readPixel(inputFile);
                tempNode->topLeft->bottomLeft = tempNode->topLeft->bottomRight = tempNode->topLeft->topLeft = tempNode->topLeft->topRight = NULL;
            } else {
                // se poate observa ca nodurile interne sunt adaugate in coada 
                tempNode->topLeft->rgbValues = NULL;
                enqueue(myQueue, tempNode->topLeft);
            }

            fread(&nodeType, 1, 1, inputFile);
            if (nodeType) {
                tempNode->topRight->rgbValues = malloc(sizeof(pixel));
                *(tempNode->topRight->rgbValues) = readPixel(inputFile);
                tempNode->topRight->bottomLeft = tempNode->topRight->bottomRight = tempNode->topRight->topLeft = tempNode->topRight->topRight = NULL;
            } else {
                tempNode->topRight->rgbValues = NULL;
                enqueue(myQueue, tempNode->topRight);
            }

            fread(&nodeType, 1, 1, inputFile);
            if (nodeType) {
                tempNode->bottomRight->rgbValues = malloc(sizeof(pixel));
                *(tempNode->bottomRight->rgbValues) = readPixel(inputFile);
                tempNode->bottomRight->bottomLeft = tempNode->bottomRight->bottomRight = tempNode->bottomRight->topLeft = tempNode->bottomRight->topRight = NULL;
            } else {
                tempNode->bottomRight->rgbValues = NULL;
                enqueue(myQueue, tempNode->bottomRight);
            }

            fread(&nodeType, 1, 1, inputFile);
            if (nodeType) {
                tempNode->bottomLeft->rgbValues = malloc(sizeof(pixel));
                *(tempNode->bottomLeft->rgbValues) = readPixel(inputFile);
                tempNode->bottomLeft->bottomLeft = tempNode->bottomLeft->bottomRight = tempNode->bottomLeft->topLeft = tempNode->bottomLeft->topRight = NULL;
            } else {
                tempNode->bottomLeft->rgbValues = NULL;
                enqueue(myQueue, tempNode->bottomLeft);
            }
        }
    }
    
    pixel **myImage = mallocPixelMatrix(size);

    initialize(myImage, size, myTree, 0, 0);
    unsigned int value = 255;

    fprintf(outputFile, "P6\n");
    fprintf(outputFile, "%u %u\n", size, size);
    fprintf(outputFile, "%u\n", value);

    // se introduc datele din matricea de pixeli in fisierul output
    int i, j;
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            writePixel(&myImage[i][j], outputFile);
        }
    }

    destroyQueue(myQueue);
    freeTree(myTree);
    freePixelMatrix(myImage, size);
}

int main(int argc, char *argv[])
{
    if (strcmp(argv[1], "-d") == 0) {
        // cerinta 3 
        FILE *inputFile = fopen(argv[2], "rb");
        FILE *outputFile = fopen(argv[3], "wb");

        decomprimare(inputFile, outputFile);

        fclose(inputFile);
        fclose(outputFile);
        
    } else if (strcmp(argv[1], "-c1") == 0) {
        // cerinta 1
        unsigned int size;

        FILE *inputFile = fopen(argv[3], "rb");
        FILE *outputFile = fopen(argv[4], "w");
        pixel** myImage = readDataFromPPM(inputFile, &size);

        unsigned int depth = getLog2(size); // adancimea maxima a arborelui

        unsigned int maxDepth = depth;
        unsigned int numOfLeaves = 0; // numarul de frunze din arbore
        unsigned int maxBlockSize = 1; // dimensiunea celui mai mare bloc al carui mean <= factor
        // valorile celor 3 variabile de mai sus se vor modifica pe parcursul functiei de construire a arborelui

        imgNode *myTree = recursiveBuild(myImage, size, 0, 0, atoi(argv[2]), depth, &maxDepth, &numOfLeaves, &maxBlockSize);
        maxDepth = depth - maxDepth + 1; // adancimea arborelui
        
        fprintf(outputFile, "%u\n", maxDepth);
        fprintf(outputFile, "%u\n", numOfLeaves);
        fprintf(outputFile, "%u\n", maxBlockSize);

        freeTree(myTree);
        freePixelMatrix(myImage, size);

        fclose(inputFile);
        fclose(outputFile);
    } else {
        // cerinta 2
        unsigned int size;

        FILE *inputFile = fopen(argv[3], "rb");
        FILE *outputFile = fopen(argv[4], "wb");
        pixel** myImage = readDataFromPPM(inputFile, &size);

        unsigned int depth = getLog2(size);
        unsigned int maxDepth = depth;
        unsigned int numOfLeaves = 0; 
        unsigned int maxBlockSize = 1;

        imgNode *myTree = recursiveBuild(myImage, size, 0, 0, atoi(argv[2]), depth, &maxDepth, &numOfLeaves, &maxBlockSize);
        maxDepth = depth - maxDepth + 1;
        
        Queue *myQueue = initQueue(); 
        levelOrder(myTree, myQueue, outputFile, size);

        destroyQueue(myQueue);
        freeTree(myTree);
        freePixelMatrix(myImage, size);
        
        fclose(inputFile);
        fclose(outputFile);
    }

    return 0;
}