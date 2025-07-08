// Vertex Shader Customizado
// Atributos de entrada
attribute vec3 position;
attribute vec2 uv;

// Uniforms
uniform mat4 projectionMatrix;
uniform mat4 modelViewMatrix;
uniform float time;

// Variáveis para passar ao fragment shader
varying vec2 vUv;
varying vec3 vPosition;

void main() {
    vUv = uv;
    vPosition = position;
    
    // Adiciona uma pequena ondulação baseada no tempo
    vec3 newPosition = position;
    newPosition.z += sin(position.x * 10.0 + time) * 0.1;
    newPosition.z += cos(position.y * 10.0 + time) * 0.1;
    
    gl_Position = projectionMatrix * modelViewMatrix * vec4(newPosition, 1.0);
} 