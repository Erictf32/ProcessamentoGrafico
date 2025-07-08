// Fragment Shader Customizado
precision mediump float;

// Uniforms
uniform float time;
uniform vec3 color;

// Variáveis vindas do vertex shader
varying vec2 vUv;
varying vec3 vPosition;

void main() {
    // Cria um padrão de cores baseado na posição e tempo
    float r = sin(vPosition.x * 5.0 + time) * 0.5 + 0.5;
    float g = cos(vPosition.y * 5.0 + time) * 0.5 + 0.5;
    float b = sin(vPosition.z * 5.0 + time * 2.0) * 0.5 + 0.5;
    
    // Adiciona efeito baseado nas coordenadas UV
    float pattern = sin(vUv.x * 20.0 + time) * cos(vUv.y * 20.0 + time);
    
    // Combina as cores
    vec3 finalColor = vec3(r, g, b) * (0.7 + pattern * 0.3);
    
    gl_FragColor = vec4(finalColor, 1.0);
} 