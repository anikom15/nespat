/* NES pattern generator
 * Copyright 2020 Ivan M. Martinez
 * MIT License */

#include <math.h>
#include <stdio.h>

#define NTSC_COLORBURST_FREQUENCY (NES_CLOCK_FREQUENCY / 6.0)
#define TIME_STEP (1.0 / (2.0 * NES_CLOCK_FREQUENCY))
#define PI 3.14159265358979323846

#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a, b) (((a) < (b)) ? (a) : (b))
#endif

const double NTSC_SYNC_LEVEL = -40.0 / 140.0;
const double NTSC_CHROMA_VPP = 40.0 / 140.0;
const double NTSC_WHITE_LEVEL = 100.0 / 140.0;
const double NTSC_BLACK_LEVEL = 7.5 / 140.0;

const double NES_CLOCK_FREQUENCY = 21.47727273e6;
const double NES_SYNC_LEVEL = 0.048;
const double NES_BLACK_LEVEL = 0.312;
const double NES_COLORBURST_LEVEL_LOW = 0.148;
const double NES_COLORBURST_LEVEL_HIGH = 0.524;
const double NES_LEVELS[2][4] = {
	{ 0.616, 0.840, 1.100, 1.100 },
	{ 0.228, 0.312, 0.552, 0.880 }
};

double g_buf_chroma_pb[41];
double g_buf_chroma_pr[41];
double g_buf_chroma_g[41];
double g_buf_notch[41];
double g_buf_luma[32];

double nes_output[524][3410];

double gamma(double x);
double ungamma(double x);
void init_nes_output(void);
double iir_filter(double x, const int N, double *, double *, const double *, const double *);
double fir_filter(double x, const int N, double *, const double *);
double chroma_filter(double x, double *buf);
double notch_filter(double x, double *buf);
double luma_filter(double x, double *buf);

int main(int argc, char *argv[])
{
	int i, j;
	double y_gain;
	double c_gain;
	double t = -7.5 * TIME_STEP;

	init_nes_output();
	y_gain = NTSC_SYNC_LEVEL / (NES_SYNC_LEVEL - NES_BLACK_LEVEL);
	c_gain = NTSC_CHROMA_VPP / (NES_COLORBURST_LEVEL_HIGH - NES_COLORBURST_LEVEL_LOW);
	printf("P3\n3410 524\n255\n");
	for (i = 0; i < 524; ++i) {
		for (j = 0; j < 3410; ++j) {
			double level;

			
			level = (nes_output[i][j] - NTSC_BLACK_LEVEL - NES_BLACK_LEVEL)
				/ (NTSC_WHITE_LEVEL - NTSC_BLACK_LEVEL);
			{
				double r, g, b;
				double y, pr, pb;

				/* m -> notch filter -> Y */
				level = luma_filter(level, g_buf_luma);
				y = notch_filter(level, g_buf_notch);

#if 0
				/* m -> mixer -> LP -> Pb */
				pb = 2.03206 * -sin(2.0 * PI * NTSC_COLORBURST_FREQUENCY * t) * level;
				pb = chroma_filter(pb, g_buf_chroma_pb);

				/* m -> mixer -> LP -> Pr */
				pr = 1.13988 * -cos(2.0 * PI * NTSC_COLORBURST_FREQUENCY * t) * level;
				pr = chroma_filter(pr, g_buf_chroma_pr);

				r = y + pr;
				g = y - 0.1942 * pb - 0.5094 * pb;
				b = y + pb;
#endif
				r = 1.629 * -sin(2.0 * PI * NTSC_COLORBURST_FREQUENCY * t + 87.649 * PI / 180.0) * level;
				r = y + chroma_filter(r, g_buf_chroma_pr);
				b = 2.451 * -sin(2.0 * PI * NTSC_COLORBURST_FREQUENCY * t + 3.54 * PI / 180.0) * level;
				b = y + chroma_filter(b, g_buf_chroma_pb);
				g = 0.622 * -sin(2.0 * PI * NTSC_COLORBURST_FREQUENCY * t + 243.717 * PI / 180.0) * level;
				g = y + chroma_filter(g, g_buf_chroma_g);
				{
					double X[3], Y[3], Z[3];
					int i;

					r = ungamma(r);
					g = ungamma(g);
					b = ungamma(b);

					Y[0] = 0.2311 * r;
					Y[1] = 0.6907 * g;
					Y[2] = 0.0782 * b;
					X[0] = 0.670 * (Y[0] / 0.330);
					X[1] = 0.265 * (Y[1] / 0.585);
					X[2] = 0.155 * (Y[2] / 0.061);
					Z[0] = (1.0 - 0.670 - 0.330) * (Y[0] / 0.330);
					Z[1] = (1.0 - 0.265 - 0.585) * (Y[1] / 0.585);
					Z[2] = (1.0 - 0.155 - 0.061) * (Y[2] / 0.061);
				
					r = 0.0;
					g = 0.0;
					b = 0.0;	
					for (i = 0; i < 3; ++i) {
						r += 3.2406 * X[i] - 1.5372 * Y[i] - 0.4986 * Z[i];
						g += -0.9689 * X[i] + 1.8758 * Y[i] + 0.0415 * Z[i];
						b += 0.0557 * X[i] - 0.2040 * Y[i] + 1.0570 * Z[i];
					}
					r = gamma(r);
					g = gamma(g);
					b = gamma(b);
				}

				r *= 255.0;
				g *= 255.0;
				b *= 255.0;

				r = min(max(r, 0.0), 255.0);
				g = min(max(g, 0.0), 255.0);
				b = min(max(b, 0.0), 255.0);

				printf("%1.f %1.f %1.f ", r, g, b);
			}
			t += TIME_STEP;
		}
		putchar('\n');
	}
	printf("\n");
	return 0;
}

double gamma(double x)
{
	if (x < 0.018)
		return 4.500 * x;
	else
		return pow(1.099 * x, 0.4500) - 0.099;
}

double ungamma(double x)
{
	if (x < 0.0812)
		return x / 4.500;
	else
		return pow((x + 0.099) / 1.099, 1.0 / 0.4500);
}

void init_nes_output(void)
{
	int i;
	int j;
	int n;
	int px, py;

	n = 0;
	px = 0;
	py = 0;
	for (i = 0; i < 524; ++i) {
		if (i < 262)
			py = i % 262 / 60;
		else
			py = (i - 262) % 262 / 60;
		for (j = 0; j < 3410; ++j) {
			px = (j - 40) / 160;
			if (i < 240 || (i >= 262 && i < 502)) {
				if (j < 40)
					nes_output[i][j] = NES_BLACK_LEVEL;
				else if (j < 2600) {
					int x;

					if (px == 0) {
						nes_output[i][j] = NES_LEVELS[0][py];
					} else if (px == 1) {
						x = (n + 6) % 12;

						nes_output[i][j] = NES_LEVELS[x < 6][py];
					} else if (px == 2) {
						x = (n + 7) % 12;

						nes_output[i][j] = NES_LEVELS[x < 6][py];
					} else if (px == 3) {
						x = (n + 8) % 12;

						nes_output[i][j] = NES_LEVELS[x < 6][py];
					} else if (px == 4) {
						x = (n + 9) % 12;

						nes_output[i][j] = NES_LEVELS[x < 6][py];
					} else if (px == 5) {
						x = (n + 10) % 12;

						nes_output[i][j] = NES_LEVELS[x < 6][py];
					} else if (px == 6) {
						x = (n + 11) % 12;

						nes_output[i][j] = NES_LEVELS[x < 6][py];
					} else if (px == 7) {
						x = n % 12;

						nes_output[i][j] = NES_LEVELS[x < 6][py];
					} else if (px == 8) {
						x = (n + 1) % 12;

						nes_output[i][j] = NES_LEVELS[x < 6][py];
					} else if (px == 9) {
						x = (n + 2) % 12;

						nes_output[i][j] = NES_LEVELS[x < 6][py];
					} else if (px == 10) {
						x = (n + 3) % 12;

						nes_output[i][j] = NES_LEVELS[x < 6][py];
					} else if (px == 11) {
						x = (n + 4) % 12;

						nes_output[i][j] = NES_LEVELS[x < 6][py];
					} else if (px == 12) {
						x = (n + 5) % 12;

						nes_output[i][j] = NES_LEVELS[x < 6][py];
					} else if (px == 13) {
						nes_output[i][j] = NES_LEVELS[1][py];
					} else {
						nes_output[i][j] = NES_LEVELS[1][1];
					}
				}
				else if (j < 2710)
					nes_output[i][j] = NES_BLACK_LEVEL;
				else if (j < 2800)
					nes_output[i][j] = NES_BLACK_LEVEL;
				else if (j < 3050)
					nes_output[i][j] = NES_SYNC_LEVEL;
				else if (j < 3090)
					nes_output[i][j] = NES_BLACK_LEVEL;
				else if (j < 3240) {
					int x;

					x = (n + 1) % 12;
					if (x < 6)
						nes_output[i][j] = NES_COLORBURST_LEVEL_LOW;
					else
						nes_output[i][j] = NES_COLORBURST_LEVEL_HIGH;
				} else if (j < 3290)
					nes_output[i][j] = NES_BLACK_LEVEL;
				else if (j < 3300)
					nes_output[i][j] = NES_BLACK_LEVEL;
				else
					nes_output[i][j] = NES_BLACK_LEVEL;
			} else if (i < 242 || (i >= 262 && i < 504)) {
				if (j < 2710)
					nes_output[i][j] = NES_BLACK_LEVEL;
				else if (j < 2800)
					nes_output[i][j] = NES_BLACK_LEVEL;
				else if (j < 3050)
					nes_output[i][j] = NES_SYNC_LEVEL;
				else if (j < 3090)
					nes_output[i][j] = NES_BLACK_LEVEL;
				else if (j < 3240) {
					int x;

					x = (n + 1) % 12;
					if (x < 6)
						nes_output[i][j] = NES_COLORBURST_LEVEL_LOW;
					else
						nes_output[i][j] = NES_COLORBURST_LEVEL_HIGH;
				} else if (j < 3290)
					nes_output[i][j] = NES_BLACK_LEVEL;
				else if (j < 3300)
					nes_output[i][j] = NES_BLACK_LEVEL;
				else
					nes_output[i][j] = NES_BLACK_LEVEL;
			} else if (i < 245 || (i >= 262 && i < 507)) {
				if (j < 2800)
					nes_output[i][j] = NES_BLACK_LEVEL;
				else if (j < 3050)
					nes_output[i][j] = NES_SYNC_LEVEL;
				else if (j < 3090)
					nes_output[i][j] = NES_BLACK_LEVEL;
				else if (j < 3240) {
					int x;

					x = (n + 1) % 12;
					if (x < 6)
						nes_output[i][j] = NES_COLORBURST_LEVEL_LOW;
					else
						nes_output[i][j] = NES_COLORBURST_LEVEL_HIGH;
				} else
					nes_output[i][j] = NES_BLACK_LEVEL;
			} else if (i < 248 || (i >= 262 && i < 510)) {
				if (j < 2570)
					nes_output[i][j] = NES_SYNC_LEVEL;
				if (j < 2800)
					nes_output[i][j] = NES_BLACK_LEVEL;
				else
					nes_output[i][j] = NES_SYNC_LEVEL;
			} else {
				if (j < 2800)
					nes_output[i][j] = NES_BLACK_LEVEL;
				else if (j < 3050)
					nes_output[i][j] = NES_SYNC_LEVEL;
				else if (j < 3090)
					nes_output[i][j] = NES_BLACK_LEVEL;
				else if (j < 3230 || (j < 3240 && i >= 262)) {
					int x;

					x = (j + 1) % 12;
					if (x < 6)
						nes_output[i][j] = NES_COLORBURST_LEVEL_LOW;
					else
						nes_output[i][j] = NES_COLORBURST_LEVEL_HIGH;
				} else
					nes_output[i][j] = NES_BLACK_LEVEL;
			}
			++n;
		}
	}
}

double chroma_filter(double x, double *buf)
{
	static const double NUM[] = {
		0.001326406105113180,
		 0.001468085009972372,
		 0.001180389895384289,
		 389.8125272098749860E-6,
		-862.3325379897424910E-6,
		-0.002390225983013109,
		-0.003856155396407588,
		-0.004791213358731894,
		-0.004645231961038205,
		-0.002861632499598719,
		 0.001032426864531606,
		 0.007334282417465559,
		 0.016087969925655622,
		 0.027035218301907322,
		 0.039606498803725895,
		 0.052957539655626988,
		 0.066049478796532085,
		 0.077763451419383373,
		 0.087034221792538682,
		 0.092983632634907443,
		 0.095033928748691904,
		 0.092983632634907443,
		 0.087034221792538682,
		 0.077763451419383373,
		 0.066049478796532085,
		 0.052957539655626988,
		 0.039606498803725895,
		 0.027035218301907322,
		 0.016087969925655622,
		 0.007334282417465559,
		 0.001032426864531606,
		-0.002861632499598719,
		-0.004645231961038205,
		-0.004791213358731894,
		-0.003856155396407588,
		-0.002390225983013109,
		-862.3325379897424910E-6,
		 389.8125272098749860E-6,
		 0.001180389895384289,
		 0.001468085009972372,
		 0.001326406105113180
	};
	return fir_filter(x, 41, buf, NUM);
}

double notch_filter(double x, double *buf)
{
	static const double NUM[] = {
		0.007529953009356627,
		 0.015773884689550526,
		 0.021568623350248600,
		 0.021742185199778136,
		 0.014404307062641945,
		 31.19176557561327810E-6,
		-0.018160911965082412,
		-0.034836799395099587,
		-0.044088977197707782,
		-0.041465443203443850,
		-0.025784378660614976,
		-31.63445949100008650E-6,
		 0.029103020939245425,
		 0.053038522216498796,
		 0.063956010385293027,
		 0.057435469376382990,
		 0.034156474980993154,
		 13.39347289969251340E-6,
		-0.035484960812888905,
		-0.062069389291461288,
		 0.928095123503055230,
		-0.062069389291461288,
		-0.035484960812888905,
		 13.39347289969251340E-6,
		 0.034156474980993154,
		 0.057435469376382990,
		 0.063956010385293027,
		 0.053038522216498796,
		 0.029103020939245425,
		-31.63445949100008650E-6,
		-0.025784378660614976,
		-0.041465443203443850,
		-0.044088977197707782,
		-0.034836799395099587,
		-0.018160911965082412,
		 31.19176557561327810E-6,
		 0.014404307062641945,
		 0.021742185199778136,
		 0.021568623350248600,
		 0.015773884689550526,
		 0.007529953009356627
	};
	return fir_filter(x, 41, buf, NUM);
}

double luma_filter(double x, double *buf)
{
	static const double NUM[] = {
		126.1265321483154200E-6,
		  0.002122357290370519,
		  0.001971527752435553,
		 -0.001206240381159128,
		 -0.004438649147019474,
		 -0.002999618795644935,
		  0.004094254219637257,
		  0.010515057904754793,
		  0.006313653857581611,
		 -0.011885863540187516,
		 -0.032399112479050807,
		 -0.031234087535086726,
		  0.012276110570105990,
		  0.096948106841850068,
		  0.193422975128580421,
		  0.257792962103390644,
		  0.257792962103390644,
		  0.193422975128580421,
		  0.096948106841850068,
		  0.012276110570105990,
		 -0.031234087535086726,
		 -0.032399112479050807,
		 -0.011885863540187516,
		  0.006313653857581611,
		  0.010515057904754793,
		  0.004094254219637257,
		 -0.002999618795644935,
		 -0.004438649147019474,
		 -0.001206240381159128,
		  0.001971527752435553,
		  0.002122357290370519,
		126.1265321483154200E-6
	};
	return fir_filter(x, 32, buf, NUM);
}

double iir_filter(double x, const int N, double *ibuf, double *obuf, const double *NUM, const double *DEN)
{
	int i;
	double y = 0.0;

	for (i = N - 1; i > 0; --i) {
		obuf[i] = obuf[i-1];
		ibuf[i] = ibuf[i-1];
	}
	ibuf[0] = x;
	for (i = 0; i < N; ++i)
		y += ibuf[i] * NUM[i];
	for (i = 1; i < N; ++i)
		y -= obuf[i] * DEN[i];
	obuf[0] = DEN[0] * y;
	return obuf[0];
}

double fir_filter(double x, const int N, double *buf, const double *NUM)
{
	int i;
	double y = 0.0;
	
	for (i = N - 1; i > 0; --i)
		buf[i] = buf[i - 1];
	buf[0] = x;
	for (i = 0; i < N; ++i)
		y += buf[i] * NUM[i];
	return y;
}
