#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lv2.h>
#include "Distortion_DS1.h"
#include "OverSample.h"

/**********************************************************************************************************************************************************/

#define PLUGIN_URI "http://moddevices.com/plugins/mod-devel/DS1"
#define TAMANHO_DO_BUFFER 256

#define INPUT_GAIN 1
#define OUTPUT_GAIN 1

enum {IN, OUT_1, TONE, LEVEL, DIST, PLUGIN_PORT_COUNT};

/**********************************************************************************************************************************************************/

class Distortion
{
public:
    Distortion() {}
    ~Distortion() {}
    static LV2_Handle instantiate(const LV2_Descriptor* descriptor, double samplerate, const char* bundle_path, const LV2_Feature* const* features);
    static void activate(LV2_Handle instance);
    static void deactivate(LV2_Handle instance);
    static void connect_port(LV2_Handle instance, uint32_t port, void *data);
    static void run(LV2_Handle instance, uint32_t n_samples);
    static void cleanup(LV2_Handle instance);
    static const void* extension_data(const char* uri);
    float *in;
    float *out_1;
    float *tone;
    float *level;
    float *dist;
    
    float *u_f;
    float *y_f;
    
    float *u;
    float *u2;
    float *y;
    float *y2;
    float *v1;
    float *v2;
    float *v3;
    
    float T;
    float SampleRate;
    
    float h1u_1;
    float h1y_1;
    
    float h2u_1;
    float h2y_1;
    float h2u_2;
    float h2y_2;
    float h2u_3;
    float h2y_3;
    float h2u_4;
    float h2y_4;
    
    float h3u_1;
    float h3y_1;
    float h3u_2;
    float h3y_2;
    
    float u_1;
    float y_1;
    float v1_1;
    float v2_1;
    float v3_1;
    
    int cont;
    
    ClipClass obj;
};

/**********************************************************************************************************************************************************/

static const LV2_Descriptor Descriptor = {
    PLUGIN_URI,
    Distortion::instantiate,
    Distortion::connect_port,
    Distortion::activate,
    Distortion::run,
    Distortion::deactivate,
    Distortion::cleanup,
    Distortion::extension_data
};

/**********************************************************************************************************************************************************/

LV2_SYMBOL_EXPORT
const LV2_Descriptor* lv2_descriptor(uint32_t index)
{
    if (index == 0) return &Descriptor;
    else return NULL;
}

/**********************************************************************************************************************************************************/

LV2_Handle Distortion::instantiate(const LV2_Descriptor* descriptor, double samplerate, const char* bundle_path, const LV2_Feature* const* features)
{
    //printf("--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**sample: %f\n", samplerate );
    Distortion *plugin = new Distortion();
    
    plugin->cont = 0;
    
    plugin->u_f = (float*)malloc(2*TAMANHO_DO_BUFFER*sizeof(float));
    plugin->y_f = (float*)malloc(2*TAMANHO_DO_BUFFER*sizeof(float));
    plugin->u = (float*)malloc(2*TAMANHO_DO_BUFFER*sizeof(float));
    plugin->y = (float*)malloc(2*TAMANHO_DO_BUFFER*sizeof(float));
    plugin->u2 = (float*)malloc(8*TAMANHO_DO_BUFFER*sizeof(float));
    plugin->y2 = (float*)malloc(8*TAMANHO_DO_BUFFER*sizeof(float));
    plugin->v1 = (float*)malloc(8*TAMANHO_DO_BUFFER*sizeof(float));
    plugin->v2 = (float*)malloc(8*TAMANHO_DO_BUFFER*sizeof(float));
    plugin->v3 = (float*)malloc(8*TAMANHO_DO_BUFFER*sizeof(float));
    
    plugin->T = 1/samplerate;
    plugin->SampleRate = samplerate;
    
    plugin->h1u_1 = 0;
    plugin->h1y_1 = 0;
    
    plugin->h2u_1 = 0;
    plugin->h2y_1 = 0;
    plugin->h2u_2 = 0;
    plugin->h2y_2 = 0;
    plugin->h2u_3 = 0;
    plugin->h2y_3 = 0;
    plugin->h2u_4 = 0;
    plugin->h2y_4 = 0;
    
    plugin->h3u_1 = 0;
    plugin->h3y_1 = 0;
    plugin->h3u_2 = 0;
    plugin->h3y_2 = 0;
    
    plugin->u_1 = 0;
    plugin->y_1 = 0;
    plugin->v1_1 = 0;
    plugin->v2_1 = 0;
    plugin->v3_1 = 0;
    	
    return (LV2_Handle)plugin;
}

/**********************************************************************************************************************************************************/

void Distortion::activate(LV2_Handle instance)
{
    // TODO: include the activate function code here
}

/**********************************************************************************************************************************************************/

void Distortion::deactivate(LV2_Handle instance)
{
    // TODO: include the deactivate function code here
}

/**********************************************************************************************************************************************************/

void Distortion::connect_port(LV2_Handle instance, uint32_t port, void *data)
{
    Distortion *plugin;
    plugin = (Distortion *) instance;

    switch (port)
    {
        case IN:
            plugin->in = (float*) data;
            break;
        case OUT_1:
            plugin->out_1 = (float*) data;
            break;
        case TONE:
            plugin->tone = (float*) data;
            break;
        case LEVEL:
            plugin->level = (float*) data;
            break;
        case DIST:
            plugin->dist = (float*) data;
            break;
    }
}

/**********************************************************************************************************************************************************/

void Distortion::run(LV2_Handle instance, uint32_t n_samples)
{
    Distortion *plugin;
    plugin = (Distortion *) instance;
    
    if( (n_samples > TAMANHO_DO_BUFFER) && (plugin->cont == 0) )
	{
		plugin->u_f = (float*)realloc(plugin->u_f, 2*n_samples*sizeof(float));
		plugin->y_f = (float*)realloc(plugin->y_f, 2*n_samples*sizeof(float));
		plugin->u = (float*)realloc(plugin->u, 2*n_samples*sizeof(float));
		plugin->y = (float*)realloc(plugin->y, 2*n_samples*sizeof(float));
		plugin->u2 = (float*)realloc(plugin->u2, 8*n_samples*sizeof(float));
		plugin->y2 = (float*)realloc(plugin->y2, 8*n_samples*sizeof(float));
		plugin->v1 = (float*)realloc(plugin->v1, 8*n_samples*sizeof(float));
		plugin->v2 = (float*)realloc(plugin->v2, 8*n_samples*sizeof(float));
		plugin->v3 = (float*)realloc(plugin->v3, 8*n_samples*sizeof(float));
    
		plugin->cont = 1;
    
		return;
	}
    
    float Tone;
    float Level;
    float Dist;
    
    Tone = (float)(*(plugin->tone));
    Level = (float)(*(plugin->level));
    Dist = (float)(*(plugin->dist));
    
    float T2;
    float T3;
    
    float SampleRate2;
    
    uint32_t n2;
    uint32_t n3;
    
    for (uint32_t i=1; i<=n_samples; i++)
    {
		plugin->in[i-1] = INPUT_GAIN*plugin->in[i-1]; //15dB
	}
	
	//Over 2x
	
	T2 = 0.5*plugin->T;
	SampleRate2 = 2*plugin->SampleRate;
    Over2_Float(plugin->in, plugin->u_f, &plugin->h1u_1, n_samples);
    n2 = 2*n_samples;
    
    /*****************************************************************/
    
    if (plugin->SampleRate == 48000)
    {
		Filter1_48000(plugin->u_f, plugin->y_f, n2, &plugin->h1u_1, &plugin->h1y_1 );
	}
	else
	{
		Filter1(plugin->u_f, plugin->y_f, n2, SampleRate2, &plugin->h1u_1, &plugin->h1y_1 );
	}
	
    /*****************************************************************/
   

    
    for (uint32_t i=1; i<=n2; i++)
    {
		plugin->u_f[i-1] = plugin->y_f[i-1]; 
	}
	
	/*****************************************************************/
	
	if (plugin->SampleRate == 48000)
    {
		Filter2_48000(plugin->u_f, plugin->y_f, n2, &plugin->h2u_1, &plugin->h2y_1, &plugin->h2u_2, &plugin->h2y_2, &plugin->h2u_3, &plugin->h2y_3, &plugin->h2u_4, &plugin->h2y_4 );
	}
	else
	{
		Filter2(plugin->u_f, plugin->y_f, n2, SampleRate2, &plugin->h2u_1, &plugin->h2y_1, &plugin->h2u_2, &plugin->h2y_2, &plugin->h2u_3, &plugin->h2y_3, &plugin->h2u_4, &plugin->h2y_4 );
	}
	
	
   /*****************************************************************/
    
   
    
    for (uint32_t i=1; i<=n2; i++)
    {
		plugin->u_f[i-1] = plugin->y_f[i-1];
	}
	
    if (plugin->SampleRate == 48000)
    {
		FilterGain_48000(plugin->u_f, plugin->y_f, n2, Dist, &plugin->h3u_1, &plugin->h3y_1, &plugin->h3u_2, &plugin->h3y_2 );
	}
	else
	{
		FilterGain(plugin->u_f, plugin->y_f, n2, Dist, SampleRate2, &plugin->h3u_1, &plugin->h3y_1, &plugin->h3u_2, &plugin->h3y_2 );
	}
	 
		
	//Over 4x
	
	T3 = 0.25*T2;
    Over4_Float(plugin->y_f, plugin->u2, &plugin->u_1, n2);
    n3 = 4*n2;



		DS1_Clip_Tone(plugin->u2, plugin->y2, plugin->v1, plugin->v2, plugin->v3, n3, T3, &plugin->u_1, &plugin->y_1, &plugin->v1_1, &plugin->v2_1, &plugin->v3_1, Tone, Level, &plugin->obj);
    
    /*****************************************************************/
    
    Down8_Float(plugin->out_1, plugin->y2, n_samples);
	
	 for (uint32_t i=1; i<=n_samples; i++)
    {
		plugin->out_1[i-1] = plugin->out_1[i-1]*OUTPUT_GAIN; //-15dB
	}
    
}

/**********************************************************************************************************************************************************/

void Distortion::cleanup(LV2_Handle instance)
{
	Distortion *plugin;
	plugin = (Distortion *) instance;
	
	free(plugin->u_f);
	free(plugin->y_f);
	free(plugin->u);
	free(plugin->y);
	free(plugin->u2);
	free(plugin->y2);
	free(plugin->v1);
	free(plugin->v2);
	free(plugin->v3);
	
	
	
    delete ((Distortion *) instance);
}

/**********************************************************************************************************************************************************/

const void* Distortion::extension_data(const char* uri)
{
    return NULL;
}
