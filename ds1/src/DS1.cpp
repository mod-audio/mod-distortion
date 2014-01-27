#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lv2.h>
#include "Distortion_DS1.h"
#include "ModFilters.h"
#include "OverSample.h"

/**********************************************************************************************************************************************************/

#define PLUGIN_URI "http://portalmod.com/plugins/mod-devel/DS1"
#define TAMANHO_DO_BUFFER 2048
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
    
    double *u;
    double *u2;
    double *y;
    double *y2;
    double *v1;
    double *v2;
    double *v3;
    
    double T;
    
    double h1u_1;
    double h1y_1;
    
    double h2u_1;
    double h2y_1;
    double h2u_2;
    double h2y_2;
    double h2u_3;
    double h2y_3;
    double h2u_4;
    double h2y_4;
    
    double h3u_1;
    double h3y_1;
    double h3u_2;
    double h3y_2;
    
    double u_1;
    double y_1;
    double v1_1;
    double v2_1;
    double v3_1;
    
    double h4u_1;
    double h4y_1;
    double h4u_2;
    double h4y_2;
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
    printf("--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**sample: %f\n", samplerate );
    Distortion *plugin = new Distortion();
    
    plugin->u = (double*)malloc(TAMANHO_DO_BUFFER*sizeof(double));
    plugin->y = (double*)malloc(TAMANHO_DO_BUFFER*sizeof(double));
    plugin->u2 = (double*)malloc(TAMANHO_DO_BUFFER*sizeof(double));
    plugin->y2 = (double*)malloc(TAMANHO_DO_BUFFER*sizeof(double));
    plugin->v1 = (double*)malloc(TAMANHO_DO_BUFFER*sizeof(double));
    plugin->v2 = (double*)malloc(TAMANHO_DO_BUFFER*sizeof(double));
    plugin->v3 = (double*)malloc(TAMANHO_DO_BUFFER*sizeof(double));
    
    plugin->T = 1/samplerate;
    
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
    
    plugin->h4u_1 = 0;
    plugin->h4y_1 = 0;
    plugin->h4u_2 = 0;
    plugin->h4y_2 = 0;
    
    	
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
    
    double Tone;
    double Level;
    double Dist;
    
    Tone = (float)(*(plugin->tone));
    Level = (float)(*(plugin->level));
    Dist = (float)(*(plugin->dist));
    
    double T2;
    double T3;
    uint32_t n2;
    uint32_t n3;
    
    for (uint32_t i=1; i<=n_samples; i++)
    {
		plugin->in[i-1] = 5.6234*plugin->in[i-1]; //15dB
	}
    
    switch (2)
    {
        case 1:
            T2 = plugin->T;
            Over1(plugin->in, plugin->u, &plugin->h1u_1, n_samples);
            n2 = n_samples;
            break;
        case 2:
            T2 = 0.5*plugin->T;
            Over2(plugin->in, plugin->u, &plugin->h1u_1, n_samples);
            n2 = 2*n_samples;
            break;
        case 3:
            T2 = 0.25*plugin->T;
            Over4(plugin->in, plugin->u, &plugin->h1u_1, n_samples);
            n2 = 4*n_samples;
            break;
        case 4:
            T2 = 0.125*plugin->T;
            Over8(plugin->in, plugin->u, &plugin->h1u_1, n_samples);
            n2 = 8*n_samples;
            break;
        default:
            T2 = plugin->T;
            Over1(plugin->in, plugin->u, &plugin->h1u_1, n_samples);
            n2 = n_samples;
    }
    
    /*****************************************************************/
    
    Filter1(plugin->u, plugin->y, n2, T2, &plugin->h1u_1, &plugin->h1y_1 );
    
    /*****************************************************************/
   

    
    for (uint32_t i=1; i<=n2; i++)
    {
		plugin->u[i-1] = plugin->y[i-1]; 
	}
	
	/*****************************************************************/
	
	Filter2(plugin->u, plugin->y, n2, T2, &plugin->h2u_1, &plugin->h2y_1, &plugin->h2u_2, &plugin->h2y_2, &plugin->h2u_3, &plugin->h2y_3, &plugin->h2u_4, &plugin->h2y_4 );
    
   /*****************************************************************/
    
   
    
    for (uint32_t i=1; i<=n2; i++)
    {
		plugin->u[i-1] = plugin->y[i-1]; 
	}
	
	
    
    FilterGain(plugin->u, plugin->y, n2, Dist, T2, &plugin->h3u_1, &plugin->h3y_1, &plugin->h3u_2, &plugin->h3y_2 );
    
   
	/*
	for (uint32_t i=1; i<=n2; i++)
    {
		plugin->u[i-1] = plugin->y[i-1]; 
	}
	*/
	
	switch (3)
    {
        case 1:
            T3 = T2;
            Over1_Double(plugin->y, plugin->u2, &plugin->h1u_1, n2);
            n3 = n2;
            break;
        case 2:
            T3 = 0.5*T2;
            Over2_Double(plugin->y, plugin->u2, &plugin->h1u_1, n2);
            n3 = 2*n2;
            break;
        case 3:
            T3 = 0.25*T2;
            Over4_Double(plugin->y, plugin->u2, &plugin->h1u_1, n2);
            n3 = 4*n2;
            break;
        case 4:
            T3 = 0.125*T2;
            Over8_Double(plugin->y, plugin->u2, &plugin->h1u_1, n2);
            n3 = 8*n2;
            break;
        default:
            T3 = T2;
            Over1_Double(plugin->y, plugin->u2, &plugin->h1u_1, n2);
            n3 = n2;
    }
    
    DS1_Clip_Tone(plugin->u2, plugin->y2, plugin->v1, plugin->v2, plugin->v3, n3, T3, &plugin->u_1, &plugin->y_1, &plugin->v1_1, &plugin->v2_1, &plugin->v3_1, Tone, Level);

    
    /*****************************************************************/
    
    switch (4)
    {
        case 1:
            Down1(plugin->out_1, plugin->y2, n_samples);
            break;
        case 2:
            Down2(plugin->out_1, plugin->y2, n_samples);
            break;
        case 3:
            Down4(plugin->out_1, plugin->y2, n_samples);
            break;
        case 4:
            Down8(plugin->out_1, plugin->y2, n_samples);
            break;
        default:
            Down1(plugin->out_1, plugin->y2, n_samples);
    }
	
    
}

/**********************************************************************************************************************************************************/

void Distortion::cleanup(LV2_Handle instance)
{
	Distortion *plugin;
	plugin = (Distortion *) instance;
	
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
