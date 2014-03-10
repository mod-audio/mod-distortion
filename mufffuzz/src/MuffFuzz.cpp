#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lv2.h>
#include "Distortion_MuffFuzz.h"
#include "ModFilters.h"
#include "OverSample.h"

/**********************************************************************************************************************************************************/

#define PLUGIN_URI "http://portalmod.com/plugins/mod-devel/MuffFuzz"
#define TAMANHO_DO_BUFFER 256
enum {IN, OUT_1, LEVEL, PLUGIN_PORT_COUNT};

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
    float *level;
    
    double *u;
    double *y;
    
    double T;
    double SampleRate;
    
    double h1u_1;
    double h1y_1;
    
    double h2u_1;
    double h2y_1;
     
    int cont;
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
    Distortion *plugin = new Distortion();
    
    plugin->cont = 0;
    
    plugin->u = (double*)malloc(2*TAMANHO_DO_BUFFER*sizeof(double));
    plugin->y = (double*)malloc(2*TAMANHO_DO_BUFFER*sizeof(double));
    
    plugin->T = 1/samplerate;
    plugin->SampleRate = samplerate;
    
    plugin->h1u_1 = 0;
    plugin->h1y_1 = 0;
    
    plugin->h2u_1 = 0;
    plugin->h2y_1 = 0;
    
    	
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
        case LEVEL:
            plugin->level = (float*) data;
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
		plugin->u = (double*)realloc(plugin->u, 2*n_samples*sizeof(double));
		plugin->y = (double*)realloc(plugin->y, 2*n_samples*sizeof(double));
    
		plugin->cont = 1;
    
		return;
	}
    
    double Level;
    
    Level = (float)(*(plugin->level));
    
    double T2;
    uint32_t n2;
    
    for (uint32_t i=1; i<=n_samples; i++)
    {
		plugin->in[i-1] = 5.6234*plugin->in[i-1]; //15dB
	}
	
	//Over 2x
	
	T2 = 0.5*plugin->T;
    Over2(plugin->in, plugin->u, &plugin->h1u_1, n_samples);
    n2 = 2*n_samples;
    
    /*****************************************************************/

	Filter1(plugin->u, plugin->y, n2, T2, &plugin->h1u_1, &plugin->h1y_1 );
	
    /*****************************************************************/
    
    for (uint32_t i=1; i<=n2; i++)
    {
		plugin->u[i-1] = plugin->y[i-1];
	}
	
	/*****************************************************************/
	
	Clip(plugin->u, plugin->y, n2);
	
   /*****************************************************************/
       
    for (uint32_t i=1; i<=n2; i++)
    {
		plugin->u[i-1] = plugin->y[i-1]; 
	}
	
	/*****************************************************************/
	
	Filter2(plugin->u, plugin->y, n2, Level, T2, &plugin->h2u_1, &plugin->h2y_1);
		
	 /*****************************************************************/
	 
	 Down2(plugin->out_1, plugin->y, n_samples);
	
	 for (uint32_t i=1; i<=n_samples; i++)
    {
		plugin->out_1[i-1] = plugin->out_1[i-1]/5.6234; //-15dB
	}
    
}

/**********************************************************************************************************************************************************/

void Distortion::cleanup(LV2_Handle instance)
{
	Distortion *plugin;
	plugin = (Distortion *) instance;
	
	free(plugin->u);
	free(plugin->y);
	
	
    delete ((Distortion *) instance);
}

/**********************************************************************************************************************************************************/

const void* Distortion::extension_data(const char* uri)
{
    return NULL;
}
